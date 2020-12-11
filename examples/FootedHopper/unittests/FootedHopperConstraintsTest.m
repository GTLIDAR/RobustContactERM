classdef FootedHopperConstraintsTest < matlab.unittest.TestCase
    %FOOTEDHOPPERCONSTRAINTSTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Optimizer object
        optimizer;
        % Test System Properties
        plant = FootedHopper();
        q;
        dq = [0.5, -0.1, 1, -2, -0.5]';
        l = [10, 0, -1, 0, 2, 1.4, 0.8, 0.2]';
        % Finite Differencing Settings
        h = 1e-8;
        tol = 1e-6;
    end
    
    methods (TestMethodSetup)
        function createOptimizer(testCase)
            % First change the foot length and inertia
            testCase.plant.lengths(3) = 1/3;
            testCase.plant.inertias(3) = 1/12 * testCase.plant.masses(3) * testCase.plant.lengths(3)^2;
            % Create the optimizer
            options = struct('frictionVariance',0.01,'heightVariance',0.1,'ermFrictionBias',0.001);
            testCase.optimizer = RobustContactImplicitTrajectoryOptimizer(testCase.plant, 101,1,options);
            % Create some seed values
            x0 = [0, 1.5, (1-testCase.plant.heelFraction)*testCase.plant.lengths(3),0]';
            testCase.q = testCase.plant.inverseKinematics(x0(1:2), x0(3:4),-30)';
        end
    end
    
    methods (Test)
        function checkIndexValues(testCase)
           % Check the normal index values
           testCase.assertEqual(testCase.optimizer.normal_inds, [1,5],'Normal force indices are incorrect');
            % Check tangential force indices
            testCase.assertEqual(testCase.optimizer.tangent_inds, [2,3,6,7], 'Tangential force indices are incorrect');
            % Check the slack variable indices
            testCase.assertEqual(testCase.optimizer.gamma_inds, [4,8], 'Sliding velocity indices are incorrect');
        end
        function normalDistanceTest(testCase)
            % Test Parameters
            x = [testCase.q; testCase.dq; testCase.l([1,5])];
            % Get the values from the normal distance 
            [f, df] = testCase.optimizer.normalDistanceConstraint(x);
            % Check the size of the distance values
            testCase.assertEqual(size(f), [2,1], 'Incorrect number of distances');
            % Check the values of the distances
            p = testCase.plant.kinematics(testCase.q);
            phi = p(2,:)';
            testCase.verifyEqual(f, phi, 'Normal distances are not accurate');
            % Estimate the derivative using finite differences
            df_est = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceConstraint(z), x);
            % Compare the values
            testCase.verifyEqual(df, df_est, 'RelTol', testCase.tol, 'Normal distance gradients are inaccurate');
        end
        function slidingVelocityTest(testCase)
            % Test Parameters
            x = [testCase.q; testCase.dq; testCase.l([1,5,4,8,2,3,6,7])];
            % Value of the sliding constraint
            [f, df] = testCase.optimizer.slidingVelocityConstraint(x);
            % Check the size of the constraint
            testCase.assertEqual(size(f), [4,1], 'Incorrect number of sliding velocity constraints');
            % Estimate the derivative using finite differences
            df_est = testCase.finiteDifference(@(z) testCase.optimizer.slidingVelocityConstraint(z), x);
            % Compare
            testCase.verifyEqual(df, df_est, 'RelTol', testCase.tol, 'Sliding velocity constraint gradients are inaccurate');
          
        end
        function slidingVelocityDefectTest(testCase)
           % Test parameters
           x = [testCase.q; testCase.dq];
           l_ = testCase.l;
           % Function values
           [~, df] = testCase.optimizer.slidingVelocityDefect(x, l_);
           % Derivative wrt state
           df_x = testCase.finiteDifference(@(z) testCase.optimizer.slidingVelocityDefect(z,l_), x);
           testCase.verifyEqual(df(:,1:numel(x)), df_x, 'RelTol',testCase.tol, 'Sliding Defect State Gradient inaccurate');
           % Derivative wrt force
           df_l = testCase.finiteDifference(@(z) testCase.optimizer.slidingVelocityDefect(x, z), l_);
           testCase.verifyEqual(df(:,numel(x)+1:end), df_l, 'RelTol', testCase.tol, 'Sliding Defect Force Gradient inaccurate');
           
        end
        function slidingVelocityValueTest(testCase)
           q_ = [1.5,2, pi/6, -pi/3, pi/2]';
           %Copy the plant and change some values
           model = testCase.plant;
           model.lengths(3) = 1/3;
           model.heelFraction = 1/4;
           
           % Check the value of the friction components of the contact
           % jacobian
           [~, ~, ~, ~, ~, ~, ~, ~, ~, Jt] = model.contactConstraints(q_);
           J_exp = [1, 0, sqrt(3)+1/8, sqrt(3)/2+1/8, +1/8;
                    1, 0, sqrt(3)-1/24, sqrt(3)/2-1/24, -1/24];
           testCase.verifyEqual(Jt{1}, J_exp,'RelTol',1e-12,'Positive Friction Vectors Inaccurate');
           testCase.verifyEqual(Jt{2}, -J_exp,'RelTol',1e-12, 'Negative Friction Vectors Inaccurate');
            
           % Check the value of the sliding velocity constraints
           dq_ = [2, -1, -1/10, 0.5, 0.01]';
           Jt_exp = [J_exp(1,:); -J_exp(1,:); J_exp(2,:); -J_exp(2,:)];
           gamma = -0.1;
           
           slidingVal = gamma + Jt_exp*dq_;
           % The constraint value
           testOptimizer = RobustContactImplicitTrajectoryOptimizer(model, 101,[1,1]);
           f = testOptimizer.slidingVelocityDefect([q_;dq_], [1,2,-3,gamma,1.4,4.2,-2.1,gamma]');
           
           % Check 
           testCase.verifyEqual(f, slidingVal,'RelTol',1e-12,'Sliding Velocity Constraint Inaccurate');           
        end
        function evalContactJacobianTest(testCase)
           % Test Case
           x = testCase.q;
           % values
           [~, dJ] = testCase.optimizer.evalContactJacobian(x);
           % Finite difference
           dJ_est = testCase.finiteDifference(@(z) testCase.optimizer.evalContactJacobian(z), x);
           % Compare
           testCase.verifyEqual(dJ, dJ_est, 'RelTol',testCase.tol, 'Contact Jacobian Derivative Inaccurate (possible ordering issues)?');
            
        end
        function contactJacobianTest(testCase)
           % Test Case
           x = testCase.q;
           % Value of the contact jacobian
           [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, dJn, dJt] = testCase.plant.contactConstraints(x);
           % Finite difference for the normal Jacobian
           JnFun = @(w) testCase.outputWrapper9(@(z) testCase.plant.contactConstraints(z), w);
           dJn_e = testCase.finiteDifference(@(v) JnFun(v), x);
           dJn_est = zeros(size(dJn,1),size(dJn,2));
           for n = 1:numel(x)
               t = dJn_e(:,:,n)';
              dJn_est(:,n) = t(:); 
           end
           % Finite difference for the tangential Jacobian
           JtFun = @(w) testCase.outputWrapper10(@(z) testCase.plant.contactConstraints(z), w);
           dJt_est = testCase.cellFiniteDifference(@(v) JtFun(v), x);
           
           % Compare the values
           testCase.verifyEqual(dJn, dJn_est, 'RelTol',testCase.tol, 'Normal Contact Jacobian gradient inaccurate');
           testCase.verifyEqual(dJt{1}, dJt_est{1}, 'RelTol', testCase.tol, 'Tangent Contact Jacobian 1 gradient inaccurate');
           testCase.verifyEqual(dJt{2}, dJt_est{2}, 'RelTol', testCase.tol, 'Tangent Contact Jacobian 2 gradient inaccurate');
           
        end
        function frictionConeTest(testCase)
            % Test Parameters
            x = testCase.l([1,5,2,3,6,7,4,8]);
            % Constraint values
            [f, df] = testCase.optimizer.frictionConeConstraint(x);
            % Check the size
            testCase.assertEqual(size(f), [2,1], 'Incorrect number of friction cone constraints');
            % Estimate the gradient
            df_est = testCase.finiteDifference(@(z) testCase.optimizer.frictionConeConstraint(z), x);
            % Compare
            testCase.verifyEqual(df, df_est, 'RelTol',testCase.tol, 'Friction cone gradient inaccurate');
        end
        function ermDistanceTest(testCase)
            % Test parameters
            lambdaN = testCase.l([1,5]);
            x = [testCase.q; testCase.dq];
            dt  = 0.01;
            % Distance ERM value
            [~, df] = testCase.optimizer.normalDistanceERMCost(x, lambdaN);
            % Check each of the derivatives
%             df_est = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(z, x, lambdaN), dt);
%             testCase.verifyEqual(df(:,1), df_est, 'RelTol',testCase.tol, 'Distance ERM stepsize gradient inaccurate');
            df_est = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(z, lambdaN), x);
            testCase.verifyEqual(df(:,2:numel(x)+1), df_est, 'RelTol',testCase.tol,'Distance ERM state gradient inaccurate');
            df_est = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(x, z), lambdaN);
            testCase.verifyEqual(df(:,numel(x)+2:end), df_est, 'RelTol', testCase.tol, 'Distance ERM Force gradient inaccurate');
            
        end
        function ermFrictionTest(testCase)
            % Test parameters
            lambda = testCase.l;
            % Constraint Values
            [~, df] = testCase.optimizer.frictionConeERMCost(lambda);
            % Estimate the gradient
            df_est = testCase.finiteDifference(@(z) testCase.optimizer.frictionConeERMCost(z), lambda);
            % Check the components of the gradient
            testCase.verifyEqual(df(:,[1,5]), df_est(:,[1,5]), 'RelTol',testCase.tol,'Friction ERM Normal Force Gradients inaccurate');
            testCase.verifyEqual(df(:,[4,8]), df_est(:,[4,8]), 'RelTol',testCase.tol,'Friction ERM Slack Variables Gradient inaccurate');
            testCase.verifyEqual(df(:,[2,3,6,7]), df_est(:,[2,3,6,7]),'RelTol',testCase.tol,'Friction ERM Tangent Force Gradients inaccurate');            
        end
    end
      methods   
        function df = finiteDifference(obj, fun, x)
            
            f = fun(x);
            
            df =  cell(1, numel(x));
            % Shifts
            dx = obj.h/2 * eye(numel(x));
            % Central differencing
            for n = 1:numel(x)
                f1 = fun(x + dx(:,n));
                f2 = fun(x - dx(:,n));
                df{n} = (f1 - f2)./obj.h;
            end
            % Collect the outputs together
            s = size(f);
            if s(2) == 1
                dim = 2;
            else
                dim = numel(s) + 1;
            end
            df = cat(dim, df{:});
        end
        function df = cellFiniteDifference(obj, fun, x)
            f = fun(x);           
            % Holding
            df = {zeros(numel(f{1}), numel(x)), zeros(numel(f{1}), numel(x))};
            % Shifts
            dx = obj.h/2 * eye(numel(x));
            % Central differencing
            for n = 1:numel(x)
                f1 = fun(x + dx(:,n));
                f2 = fun(x - dx(:,n));
                % Transpose
                f1{1} = f1{1}';
                f1{2} = f1{2}';
                f2{1} = f2{1}';
                f2{2} = f2{2}';
                % Unwrap
                df{1}(:,n) = (f1{1}(:) - f2{1}(:))./obj.h;
                df{2}(:,n) = (f1{2}(:) - f2{2}(:))./obj.h;
            end
        end
      end
      methods (Static)
          function fun2 = outputWrapper2(fun, x)
              [~, fun2] = fun(x);
          end
          function fun9 = outputWrapper9(fun, x)
              [~, ~, ~, ~, ~, ~, ~, ~, fun9] = fun(x);
          end
          function fun10 = outputWrapper10(fun, x)
              [~, ~, ~, ~, ~, ~, ~, ~, ~, fun10] = fun(x);
          end
      end
end

