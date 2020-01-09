classdef DistanceERMTest < matlab.unittest.TestCase
    % DISTANCEERMTEST: Unit testing framework for the Gaussian ERM Class
    % UncertainDistanceERM
    %
    %   DistanceERMTest is a suite of unit tests for evaluating the outputs
    %   of the methods of the DistanceERMClass and the GaussianERM class.
    %   DistanceERMTest uses a Block as the plant model for testing the
    %   classes, although the dynamics of the plant are not used here. An
    %   example LCP problem is given for the purposes of testing the
    %   outputs of various functions. 
    %
    
    %   Luke Drnach
    %   January 9, 2020
    
    % Verified Correct (with the example)
    %   1. All derivatives of the terrain height uncertainty mean (m_mu) and
    %       deviation (m_sigma)
    %   2. Sizes of return values for all ERM method outputs
    %   3. ERM Solution is nonnegative definite
    %
    % Verified approximately correct (numerically verified)
    %   1. ERM Cost function gradient and hessian wrt the decision
    %   variables
    %   2. Gradient and Hessian of PDF and CDF values wrt decision
    %   variables
    %   3. 
    %
    % Unverified
    %   1. Cost function mixed partial derivatives
    %   2. ERM solution gradient with respect to parameters
    %   3. PDF and CDF gradients and mixed partials with respect to
    %   parameters
    
    properties
        solver;
        P;      % LCP Problem Matrix
        w;      % LCP offset vector
        x;      % LCP solution
        dP;     % LCP matrix gradient
        dw;     % LCP offset gradient
        % Finite difference parameters
        h = 1e-8;       % Step size for finite differencing
        tol = 1e-6;     % Relative Error Tolerance for finite differencing
    end
    methods (TestMethodSetup)
        function createExample(testCase)
            % Create the plant needed to create a solver
            plant = Block();    % Use a block with friction as our example
            mu = 0;             % Test value for height error
            sig = 0.1;          % Test value for height uncertainty
            testCase.solver = UncertainDistanceERM(plant,mu,sig);
            % Example LCP parameters
            testCase.P = [1    2  0 0;      %LCP problem matrix
                          2    1  0 0;
                          0    0  1 1;
                          0.5 -1 -1 0];
            testCase.w = [-2, 2, -2, 0]';   %LCP problem vector
            testCase.x = [2, 0, 1, 1]';     %LCP Solution vector
            % Example LCP Gradients
            testCase.dP = zeros(4,4,2);
            testCase.dP(:,:,1) = [1 -1 0 0; -1 2 1 0; 0 1 0 0; 0 0 0 0];
            testCase.dP(:,:,2) = [1 0 1 0; 0 0 -1 0; 1 -1 2 0; 0 0 0 0];
            testCase.dw = [0 2; 1 1; -1 -1; 0 0]; 
        end
    end
    methods (Test)
        function testMean(testCase)
            % TESTMEAN: Test the ERMMean function
            % Get the values for the mean
            [mu, dmu_x, dmu_y, dmu_xx, dmu_xy] = testCase.solver.ermMean(testCase.x,testCase.P,testCase.w,testCase.dP,testCase.dw); 
            % Check the returned value for the mean
            testCase.verifyEqual(mu, 0,'Failed value check');
            % Check the derivative wrt x
            testCase.verifyEqual(dmu_x, [1 2 0 0],'Failed solution first derivative check');
            % Check the derivative wrt y (extra parameters)
            testCase.verifyEqual(dmu_y, [2,5],'Failed parameter first derivative check');
            % Check the second derivative wrt x
            testCase.verifyEqual(dmu_xx, zeros(1,4,4),'Failed solution second derivative check');
            % Check the mixed second derivative
            testCase.verifyEqual(dmu_xy, testCase.dP(1,:,:),'Failed mixed derivative check');
        end
        function testDeviation(testCase)
            % TESTDEVIATION: Test the ERMDeviation function           
            % Get the value of the standard deviation
            [m_sigma, dsigma_x, dsigma_y, dsigma_xx, dsigma_xy]  = testCase.solver.ermDeviation(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
            % Check the returned value for the standard deviation
            testCase.verifyEqual(m_sigma, 10,'Failed value check');
            % Check the derivative wrt x
            testCase.verifyEqual(dsigma_x, zeros(1,4),'Failed first solution derivative check');
            % Check the derivative wrt y (extra parameters)
            testCase.verifyEqual(dsigma_y, zeros(1,2),'Failed first parameter derivative check');
            % Check the second derivative wrt x
            testCase.verifyEqual(dsigma_xx, zeros(1,4,4),'Failed solution second derivative check');
            % Check the mixed second derivative
            testCase.verifyEqual(dsigma_xy, zeros(1,4,2),'Failed mixed derivative check');
        end
        % --- Tests for the GaussianERM Superclass --- %
        function testCostFunctionValues(testCase)
            % Get the function value
            [r, dr] = testCase.solver.ermCost(testCase.x, testCase.P, testCase.w);
            % Check the sizes of the values as a first-pass
            testCase.verifyEqual(size(r), [1,1], 'Residual is the wrong size');
            testCase.verifyEqual(size(dr), [1, numel(testCase.x)],'Gradient is the wrong size');
            % Check for accuracy of the values using finite difference
            fun = @(z) testCase.solver.ermCost(z, testCase.P, testCase.w);
            dr_est = testCase.finiteDifference(fun, testCase.x);
            testCase.verifyEqual(dr, dr_est, 'RelTol',testCase.tol, "Gradient doesn't match numerical gradient to desired precision");
        end
        function testDistributionValues(testCase)
            nX = numel(testCase.x);
            nY = size(testCase.dP,3);
            % Get the values of the pdf and cdf
            [pdf, cdf, dp_x, dc_x, dp_y, dc_y, dp_xx, dc_xx, dp_yx, dc_yx] = testCase.solver.evalDistribution(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
            % Check that the sizes are accurate
            testCase.verifyEqual(size(pdf), [1,1],'PDF value is the wrong size');
            testCase.verifyEqual(size(cdf), [1,1],'CDF value is the wrong size');
            testCase.verifyEqual(size(dp_x), [1,nX],'dp_x value is the wrong size');
            testCase.verifyEqual(size(dc_x),[1,nX],'dc_x value is the wrong size');
            testCase.verifyEqual(size(dp_y), [1,nY],'dp_y value is the wrong size');
            testCase.verifyEqual(size(dc_y),[1, nY],'dc_y value is the wrong size');
            testCase.verifyEqual(size(dp_xx),[1, nX, nX],'dp_xx value is the wrong size');
            testCase.verifyEqual(size(dc_xx),[1, nX, nX],'dc_xx value is the wrong size');
            testCase.verifyEqual(size(dp_yx),[1,nX, nY],'dp_yx value is the wrong size');
            testCase.verifyEqual(size(dc_yx), [1,nX, nY],'dc_yx value is the wrong size');
            % Check the gradients with respect to the solutions
            p_fun = @(z) testCase.solver.evalDistribution(z, testCase.P, testCase.w);
            c_fun = @(z) testCase.outputWrapper2(p_fun, z);
            dp_est = testCase.finiteDifference(p_fun, testCase.x);
            dc_est = testCase.finiteDifference(c_fun, testCase.x);
            testCase.verifyEqual(dp_x, dp_est, 'RelTol',testCase.tol,'pdf gradient does not match numerical gradient to desired precision');
            testCase.verifyEqual(dc_x ,dc_est, 'RelTol',testCase.tol, 'cdf gradient does not match numerical gradient to desired precision');
            % Verify the second derivatives with respect to x
            dp_fun = @(z) testCase.solver.evalDistribution(z, testCase.P, testCase.w, testCase.dP, testCase.dw);
            px_fun = @(z) testCase.outputWrapper3(dp_fun, z);
            cx_fun = @(z) testCase.outputWrapper4(dp_fun, z);
            pxx_est = testCase.finiteDifference(px_fun, testCase.x);
            cxx_est  =testCase.finiteDifference(cx_fun, testCase.x);
            testCase.verifyEqual(squeeze(dp_xx), pxx_est, 'RelTol',testCase.tol, 'pdf hessian does not match numerical hessian to desired precision');
            testCase.verifyEqual(squeeze(dc_xx), cxx_est, 'RelTol', testCase.tol, 'cdf hessian does not match numerical hessian to desired precision');
        end
        function testSecondDerivatives(testCase)
           nX = numel(testCase.x);
           nY = size(testCase.dP, 3);
           % Get the values of the second derivatives of the cost
           [g_xx, g_xy] = testCase.solver.ermSecondDerivatives(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
           % Check that the sizes are correct
           testCase.verifyEqual(size(g_xx), [nX, nX], 'Hessian g_xx is the wrong size');
           testCase.verifyEqual(size(g_xy), [nX, nY], 'Mixed derivative g_xy is the wrong size');
           % Test that the second derivative matches the finite difference
           fun = @(z) testCase.solver.ermCost(z, testCase.P, testCase.w);
           fun2 = @(z) testCase.outputWrapper2(fun, z);
           g_xx_est = testCase.finiteDifference(fun2, testCase.x);
           testCase.verifyEqual(g_xx, g_xx_est, 'RelTol',testCase.tol,'Hessian g_xx does not match numerical Hessian to the desired precision');
           
        end
        function testGradient(testCase)
           nX = numel(testCase.x);
           nY = size(testCase.dP, 3);
           % Get the value of the gradient of the solution with respect to
           % the parameters
           df = testCase.solver.gradient(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
           % Check that the size is correct
           testCase.verifyEqual(size(df), [nX, nY], 'Gradient is the wrong size');
        end
        function testSolutionSize(testCase)
            % Check that the solver returns a solution
            [f, r] = testCase.solver.solve(testCase.P, testCase.w);
            % Check that the solution is the right size
            testCase.verifyEqual(size(f), size(testCase.x));
            testCase.verifyEqual(size(r), [1,1]);
            % Check that the solution is nonnegative
            testCase.verifyTrue(all(f >= 0), 'ERM solution is negative');
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
    methods (Static)
        function f = outputWrapper2(fun, x)
            [~, f] = fun(x);
        end
        function f = outputWrapper3(fun,x)
           [~, ~, f] = fun(x); 
        end
        function f = outputWrapper4(fun, x)
           [~, ~, ~, f] = fun(x); 
        end
    end
end