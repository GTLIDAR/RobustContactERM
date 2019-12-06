classdef ContactCartTests < matlab.unittest.TestCase
    
    
    properties
        plant = ContactDrivenCart();
        h = 1e-8;
        testConfig = [2, pi/3, -3*pi/4]';
        testVel = [0.1, -0.3, 0.2]';
        testU = [-0.25, 0.01]';
        testT = 0.01;
        tolerance = 1e-6;
    end
    
    methods (Test)
        % --- TEST THE PROPERTIES OF THE CART --------- %
        function massGradientTest(testCase)
            % Test the gradient of the mass matrix
            q = testCase.testConfig;
            [~, df] = testCase.plant.massMatrix(q);
            % Now calculate the gradient numerically
            fun = @(x) testCase.plant.massMatrix(x);
            df_est = testCase.finiteDifference(fun, q);
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df_est,df, 'RelTol', testCase.tolerance);
        end
        function gravityGradientTest(testCase)
            % Test the gradient of the gravity matrix
            q = testCase.testConfig;
            [~, df] = testCase.plant.gravityMatrix(q);
            % Now calculate the gradient numerically
            fun = @(x) testCase.plant.gravityMatrix(x);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(2,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df_est,df, 'RelTol', testCase.tolerance);
        end
        function jacobianGradientTest(testCase)
            % Test the gradient of the jacobian matrix
            q = testCase.testConfig;
            [~, df] = testCase.plant.jacobian(q);
            % Now calculate the gradient numerically
            fun = @(x) testCase.plant.jacobian(x);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df_est,df, 'RelTol', testCase.tolerance);
        end
        function coriolisPositionGradientTest(testCase)
            % Test the gradient of the jacobian matrix
            q = testCase.testConfig;
            dq = testCase.testVel;
            [~, df] = testCase.plant.coriolisMatrix(q,dq);
            df = df(:,:,1:3);
            % Now calculate the gradient numerically
            fun = @(x) testCase.plant.coriolisMatrix(x,dq);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df_est,df, 'RelTol', testCase.tolerance);
        end
        function coriolisVelocityGradientTest(testCase)
            % Test the gradient of the jacobian matrix
            q = testCase.testConfig;
            dq = testCase.testVel;
            [~, df] = testCase.plant.coriolisMatrix(q,dq);
            df = df(:,:,4:end);
            % Now calculate the gradient numerically
            fun = @(x) testCase.plant.coriolisMatrix(q,x);
            df_est = testCase.finiteDifference(fun, dq);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df_est,df, 'RelTol', testCase.tolerance);
        end
        function controllerGradientTest(testCase)
            % Test the gradient of the mass matrix
            q = testCase.testConfig;
            [~, df] = testCase.plant.controllerMatrix(q);
            % Now calculate the gradient numerically
            fun = @(x) testCase.plant.controllerMatrix(x);
            df_est = testCase.finiteDifference(fun, q);
            % Calculate the absolute value of the largest difference
            df_est = cat(3,df_est{:});
            % Check that the gradients are close
            testCase.verifyEqual(df_est,df, 'RelTol', testCase.tolerance);
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

