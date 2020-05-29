classdef terrainHeightUnitTest < matlab.unittest.TestCase
    
    
    properties
        optimizer;
        inContact;
        noContact;
        % Finite difference parameters
        h = 1e-6;
        tol = 1e-8;
    end
    
    methods (TestMethodSetup)
        function createExamples(testCase)
            plant = ContactDrivenCart();
            plant.timestep = 0.01;
            plant.cartHeight = 1.5;     % Change the cart height
            
            % Specify in-contact and not-in-contact conditions
            x0 = [0,0]';      % Endpoint in contact
            xf = [0,1]';      % Endpoint not in contact
            % Calculate initial and final states
            q0 = plant.inverseKinematics(x0);
            qf = plant.inverseKinematics(xf);
            % Set the desired velocities to zero
            x0 = [q0; zeros(3,1)];
            xf = [qf; zeros(3,1)];
            
            % Store these values
            testCase.inContact = struct('x',x0,'h',plant.timestep,'f',200);
            testCase.noContact = struct('x',xf,'h',plant.timestep,'f',0);
            
            % Create an optimization problem
            options.integration_method = RobustContactImplicitTrajectoryOptimizer.BACKWARD_EULER;
            options.uncertainty_source = RobustContactImplicitTrajectoryOptimizer.DISTANCE_UNCERTAINTY;
            options.heightVariance = 0.1;
            options.distribution = RobustContactImplicitTrajectoryOptimizer.GAUSSIAN;
            options.contactCostMultiplier = 1;
            
            testCase.optimizer = RobustContactImplicitTrajectoryOptimizer(plant, 101, 1, options);
        end
        
        
    end
    
    methods (Test)
        function ermGradientContactTest(testCase)
            
            h_ = testCase.inContact.h;
            x_ = testCase.inContact.x;
            l_ = testCase.inContact.f;
            [~, df] = testCase.optimizer.normalDistanceERMCost(h_, x_, l_);
            
            % Test the gradient wrt the timestep h
            df_h = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(z, x_, l_), h_);
            testCase.verifyEqual(df(:,1),df_h,'RelTol',testCase.tol,'Step size gradient inaccurate');
            
            % Test the gradient wrt the state x
            df_x = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, z, l_), x_);
            testCase.verifyEqual(df(:,2:1+length(x_)), df_x, 'RelTol',testCase.tol,'State gradient inaccurate');
            
            % Test the gradient wrt the normal force l
            df_l = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, x_, z), l_);
            testCase.verifyEqual(df(:,length(x_)+2:end), df_l, 'RelTol',testCase.tol,'Force gradient inaccurate');
        end
        function ermGradientNoContactTest(testCase)
            
            h_ = testCase.noContact.h;
            x_ = testCase.noContact.x;
            l_ = testCase.noContact.f;
            [~, df] = testCase.optimizer.normalDistanceERMCost(h_, x_, l_);
            
            % Test the gradient wrt the timestep h
            df_h = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(z, x_, l_), h_);
            testCase.verifyEqual(df(:,1),df_h,'RelTol',testCase.tol,'Step size gradient inaccurate');
            
            % Test the gradient wrt the state x
            df_x = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, z, l_), x_);
            testCase.verifyEqual(df(:,2:1+length(x_)), df_x, 'RelTol',testCase.tol,'State gradient inaccurate');
            
            % Test the gradient wrt the normal force l
            df_l = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, x_, z), l_);
            testCase.verifyEqual(df(:,length(x_)+2:end), df_l, 'RelTol',testCase.tol,'Force gradient inaccurate');
        end
    end
    
    methods
        function df = finiteDifference(obj, fun, x)
            
            % Output size
            f = fun(x);
            df = zeros(numel(f), numel(x));
            % Shifts
            dx = obj.h./2 .* eye(numel(x));
            
            % Use a central difference
            for n = 1:numel(x)
                f1 = fun(x + dx(:,n));
                f2 = fun(x - dx(:,n));
                df(:,n) = (f1 - f2)./obj.h;
            end
        end
    end
end