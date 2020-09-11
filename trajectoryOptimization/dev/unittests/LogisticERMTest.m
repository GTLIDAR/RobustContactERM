classdef LogisticERMTest < matlab.unittest.TestCase
    %LOGISTICERMTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        logisticFun;
        h = 1e-8;
        tol = 1e-6;        
        mean;
        val;
        sigma;
    end
    
    methods (TestMethodSetup)
        function createTests(testCase)
            % Get the handle to the test function
            testCase.logisticFun = @RobustContactImplicitTrajectoryOptimizer.ermCostLogistic;
            % Set the example values
            testCase.mean = [1, 2]';
            testCase.val = [0.5, 0.1]';
            testCase.sigma = [0.1,1.5]';
            
        end
    end
    methods (Test)
        function gradientTest(testCase)
           % Get the values from the function
           [f, dfx, dfmu, dfsigma] = testCase.logisticFun(testCase.val, testCase.mean, testCase.sigma); 
           % Check the size of the output
           testCase.assertEqual(size(f), size(testCase.val));
           % Check the gradient wrt the value using central differencing
           dfx_est = testCase.finiteDifference(@(x) testCase.logisticFun(x, testCase.mean, testCase.sigma), testCase.val);
           testCase.verifyEqual(dfx, dfx_est, 'RelTol',testCase.tol,'VAL gradient inaccurate');
           
           % Check the gradient wrt the mean using central differencing 
           dfmu_est = testCase.finiteDifference(@(x) testCase.logisticFun(testCase.val, x, testCase.sigma), testCase.mean);
           testCase.verifyEqual(dfmu, dfmu_est, 'RelTol',testCase.tol, 'MU Gradient inaccurate');
           
           % Check the gradient wrt sigma using central differencing
           dfs_est = testCase.finiteDifference(@(x) testCase.logisticFun(testCase.val, testCase.mean, x), testCase.sigma);
           testCase.verifyEqual(dfsigma, dfs_est, 'RelTol',testCase.tol, 'SIGMA Gradient inaccurate');
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

