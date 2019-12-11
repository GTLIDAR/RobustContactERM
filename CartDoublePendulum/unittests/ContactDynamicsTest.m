classdef ContactDynamicsTest < matlab.unittest.TestCase
   
    properties 
        plant = ContactDrivenCart();
        h = 1e-8;
        testConfig = [2, pi/3, -3*pi/4]';
        testVel = [0.1, -0.3, 0.2]';
        testU = [-0.25, 0.01]';
        testT = 0.01;
        tolerance = 1e-4;
        
        contactConfig = [0, pi/3, -2*pi/3]';
        contactVel = [0, -0.01, 0.02]';
        contactU = [-0.01, 0.025]';
    end
    methods (Test)
        % --- TEST THE PROPERTIES OF CONTACTDYNAMICS (NO CONTACT) --- %
        function dynamicsGradientTest_NoContact(testCase)
            % Get the test values
            q = testCase.testConfig;
            t = testCase.testT;
            dq = testCase.testVel;
            u = testCase.testU;
            x = [q;dq];
            % Check that the contact force is zero
            fc = testCase.plant.contactForce(q, dq,u);
            testCase.assertEqual(fc, zeros(size(fc)),'Nonzero contact force present');
            % Get the value of the dynamics
            [~, df] = testCase.plant.dynamics(t,x,u);
            %% Gradient WRT Time
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.dynamics(w, x, u);
            df_est = testCase.finiteDifference(fun, t);
            % Calculate the absolute value of the largest difference
            df_est = df_est{:};
            % Check that the gradients are close
            testCase.verifyEqual(df(:,1),df_est, 'RelTol', testCase.tolerance,'Time Derivative Inaccurate');
            %% GRADIENT WRT STATE
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.dynamics(t, w, u);
            df_est = testCase.finiteDifference(fun, x);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,2:numel(q)+1),df_est(:,1:numel(q)), 'RelTol', testCase.tolerance,'Configuration Derivative Inaccurate');
            testCase.verifyEqual(df(:,2+numel(q):numel(x)+1),df_est(:,numel(q)+1:end), 'RelTol', testCase.tolerance,'Velocity Derivative Inaccurate');
            %% GRADIENT WRT CONTROL
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.dynamics(t, x, w);
            df_est = testCase.finiteDifference(fun, u);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,numel(x)+2:end),df_est, 'RelTol', testCase.tolerance,'Control Derivative Inaccurate');
        end
        function forceGradientTest_NoContact(testCase)
            % Get the test values
            q = testCase.testConfig;
            dq = testCase.testVel;
            u = testCase.testU;
            % Get the value of the dynamics
            [fc, df] = testCase.plant.contactForce(q, dq,u);
            % Assert that the force is nonzero
            testCase.assertGreaterThanOrEqual(fc, zeros(size(fc)),'Negative LCP solution present');
            testCase.assertEqual(fc, zeros(size(fc)),'Nonzero contact force present');
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.contactForce(w, dq, u);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,1:numel(q)),df_est, 'RelTol', testCase.tolerance,'Configuration Derivative Inaccurate');
            %% Velocity Derivative
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.contactForce(q, w, u);
            df_est = testCase.finiteDifference(fun, dq);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,1+numel(q):2*numel(q)),df_est, 'RelTol', testCase.tolerance,'Velocity Derivative Inaccurate');
            %% Control Derivative
            fun = @(w) testCase.plant.contactForce(q, dq, w);
            df_est = testCase.finiteDifference(fun, u);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,2*numel(q)+1:end),df_est, 'RelTol', testCase.tolerance,'Control Derivative Inaccurate');
            
        end
        % --- TEST THE DYNAMICS GRADIENTS WITH CONTACT ---- %
        function forceGradientTest(testCase)
            % Get the test values
            q = testCase.contactConfig;
            dq = testCase.contactVel;
            u = testCase.contactU;
            % Get the value of the dynamics
            [fc, df] = testCase.plant.contactForce(q, dq,u);
            % Assert that the force is nonzero
            testCase.assertGreaterThanOrEqual(fc, zeros(size(fc)),'Negative LCP solution present');
            testCase.assertNotEmpty(fc(fc>0),'Contact force is zero');
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.contactForce(w, dq, u);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,1:numel(q)),df_est, 'RelTol', testCase.tolerance,'Configuration Derivative Inaccurate');
            %% Velocity Derivative
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.contactForce(q, w, u);
            df_est = testCase.finiteDifference(fun, dq);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,1+numel(q):2*numel(q)), df_est,'RelTol', testCase.tolerance,'Velocity Derivative Inaccurate');
            %% Control Derivative
            fun = @(w) testCase.plant.contactForce(q, dq, w);
            df_est = testCase.finiteDifference(fun, u);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,2*numel(q)+1:end),df_est, 'RelTol', testCase.tolerance,'Control Derivative Inaccurate'); 
        end
        function dynamicsGradientTest(testCase)
            % Get the test values
            q = testCase.contactConfig;
            dq = testCase.contactVel;
            u = testCase.contactU;
            t = testCase.testT;
            x = [q;dq];
            % Check that the contact force is zero
            fc = testCase.plant.contactForce(q, dq,u);
            testCase.assertGreaterThanOrEqual(fc, zeros(size(fc)),'Negative LCP solution present');
            testCase.assertNotEmpty(fc(fc>0),'Contact force is zero');
            % Get the value of the dynamics
            [~, df] = testCase.plant.dynamics(t,x,u);
            %% Gradient WRT Time
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.dynamics(w, x, u);
            df_est = testCase.finiteDifference(fun, t);
            % Calculate the absolute value of the largest difference
            df_est = df_est{:};
            % Check that the gradients are close
            testCase.verifyEqual(df(:,1),df_est, 'RelTol', testCase.tolerance,'Time Derivative Inaccurate');
            %% GRADIENT WRT STATE
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.dynamics(t, w, u);
            df_est = testCase.finiteDifference(fun, x);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,2:numel(q)+1),df_est(:,1:numel(q)), 'RelTol', testCase.tolerance,'Configuration Derivative Inaccurate');
            testCase.verifyEqual(df(:,2+numel(q):numel(x)+1),df_est(:,numel(q)+1:end), 'RelTol', testCase.tolerance,'Velocity Derivative Inaccurate');
            %% GRADIENT WRT CONTROL
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.dynamics(t, x, w);
            df_est = testCase.finiteDifference(fun, u);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,numel(x)+2:end),df_est, 'RelTol', testCase.tolerance,'Control Derivative Inaccurate');
        end        
        function contactJacobianNormalGradientTest(testCase)
            % Get the test values
            q = testCase.contactConfig;
            % Get the value of the dynamics
            %% NORMAL COMPONENT OF THE JACOBIAN
            [~, ~, dfN, dfT] = testCase.plant.contactJacobian(q);
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.contactJacobian(w);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(dfN,df_est, 'RelTol', testCase.tolerance,'Normal Contact Jacobian Gradient Inaccurate');
            %% TANGENT COMPONENT OF THE JACOBIAN
            % Now calculate the gradient numerically
            fun2 = @(w) testCase.secondOutputWrapper(fun, w);
            df_est = testCase.finiteDifference(fun2, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(dfT,df_est, 'RelTol', testCase.tolerance,'Tangent Contact Jacobian Gradient Inaccurate');
        end
        % --- TEST THE LCP PARAMETER GRADIENTS WITH CONTACT --- %
        function lcpMatrixGradientTest(testCase)
            % Get the test values
            q = testCase.contactConfig;
            dq = testCase.contactVel;
            u = testCase.contactU;
            % Get the value of the dynamics
            [~, ~, df] = testCase.plant.getLCP(q, dq,u);
            %% Position Gradient Test
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.getLCP(w, dq, u);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,:,1:numel(q)),df_est, 'RelTol', testCase.tolerance,'Configuration Derivative Inaccurate');
            %% Velocity Gradient Test
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.getLCP(q, w, u);
            df_est = testCase.finiteDifference(fun, dq);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,:,numel(q)+1:2*numel(q)),df_est, 'RelTol', testCase.tolerance,'Velocity Derivative Inaccurate');
            %% Control Gradient Test
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.getLCP(q, dq, w);
            df_est = testCase.finiteDifference(fun, u);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,:,2*numel(q)+1:end),df_est, 'RelTol', testCase.tolerance,'Control Derivative Inaccurate');
        end
        function lcpVectorGradientTest(testCase)
            % Get the test values
            q = testCase.contactConfig;
            dq = testCase.contactVel;
            u = testCase.contactU;
            % Get the value of the dynamics
            [~, ~, ~, df] = testCase.plant.getLCP(q, dq,u);
            %% Position Gradient Test
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.getLCP(w, dq, u);
            fun2 = @(w) testCase.secondOutputWrapper(fun, w);
            df_est = testCase.finiteDifference(fun2, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,1:numel(q)),df_est, 'RelTol', testCase.tolerance,'Configuration Derivative Inaccurate');
            %% Velocity Gradient Test
            fun = @(w) testCase.plant.getLCP(q, w, u);
            fun2 = @(w) testCase.secondOutputWrapper(fun, w);
            df_est = testCase.finiteDifference(fun2, dq);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,numel(q)+1:2*numel(q)),df_est, 'RelTol', testCase.tolerance,'Velocity Derivative Inaccurate');
            %% Control Gradient Test
            % Now calculate the gradient numerically
            fun = @(w) testCase.plant.getLCP(q, dq, w);
            fun2 = @(w) testCase.secondOutputWrapper(fun, w);
            df_est = testCase.finiteDifference(fun2, u);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df(:,2*numel(q)+1:end),df_est, 'RelTol', testCase.tolerance,'Control Gradient Inaccurate');
        end
    end
    methods
        function df = finiteDifference(obj,fun, q)
            
            N = length(q);
            df = cell(1,N);
            for n = 1:N
                dq = zeros(N,1);
                dq(n) = 0.5*obj.h;
                % Use a finite central difference
                f1 = fun(q + dq);
                f2 = fun(q - dq);
                % Estimate the gradient
                df{n} = (f1 - f2)./obj.h;
            end
        end
    end
    methods (Static)
        function f = secondOutputWrapper(fun, x)
            [~, f] = fun(x);
        end
    end
end