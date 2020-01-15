classdef FrictionERMTest < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Luke Drnach
    %   January 7, 2019
    
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
            mu = 0.5;           % Test value for friction coefficient
            sig = 0.1;          % Test value for friction uncertainty
            testCase.solver = UncertainFrictionERM(plant,mu,sig);
            % Example LCP parameters
            testCase.P = [1    2  0 0;      %LCP problem matrix
                          2    1  0 1;
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
        %% --- Tests for the UncertainFrictionERM Class ----- %%
        function testMean(testCase)
            %% TESTMEAN: Test the ERMMean function
                       
            % Get the values for the mean
            [mu, dmu_x, dmu_xx, dmu_y, dmu_xy] = testCase.solver.ermMean(testCase.x,testCase.P,testCase.w,testCase.dP,testCase.dw);
            
            % Check the returned value for the mean
            testCase.verifyEqual(mu, 0,'Failed value check');
            % Check the derivative wrt x
            testCase.verifyEqual(dmu_x, [0.5, -1, -1, 0],'Failed solution first derivative check');
            % Check the derivative wrt y (extra parameters)
            testCase.verifyEqual(dmu_y, [0,0],'Failed parameter first derivative check');
            % Check the second derivative wrt x
            testCase.verifyEqual(dmu_xx, zeros(1,4,4),'Failed solution second derivative check');
            % Check the mixed second derivative
            testCase.verifyEqual(dmu_xy, zeros(1,4,2),'Failed mixed derivative check');
        end
        function testDeviation(testCase)
            %% TESTDEVIATION: Test the ERMDeviation function
            
            % Get the value of the standard deviation
            [m_sigma, dsigma_x, dsigma_xx,  dsigma_y, dsigma_xy]  = testCase.solver.ermDeviation(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
            % Check the returned value for the standard deviation
            testCase.verifyEqual(m_sigma, 0.2,'Failed value check');
            % Check the derivative wrt x
            testCase.verifyEqual(dsigma_x, [0.1, 0, 0, 0],'Failed first solution derivative check');
            % Check the derivative wrt y (extra parameters)
            testCase.verifyEqual(dsigma_y, zeros(1,2),'Failed first parameter derivative check');
            % Check the second derivative wrt x
            testCase.verifyEqual(dsigma_xx, zeros(1,4,4),'Failed solution second derivative check');
            % Check the mixed second derivative
            testCase.verifyEqual(dsigma_xy, zeros(1,4,2),'Failed mixed derivative check');
        end
        %% --- Tests for the GaussianERM Superclass --- %%
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
            [pdf, cdf, dp_x, dc_x, dp_xx, dc_xx, dp_y, dc_y, dp_yx, dc_yx] = testCase.solver.evalDistribution(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
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
            % Check that distribution values are accurate
            testCase.verifyEqual(pdf, normpdf(1, 0, 0.2), 'PDF Value inaccurate');
            testCase.verifyEqual(cdf, normcdf(1, 0, 0.2), 'CDF Value inaccurate');
            % Check that the zero derivatives are zero
            testCase.verifyEqual(dp_y, zeros(1, nY), 'dp_y is nonzero');
            testCase.verifyEqual(dp_yx, zeros(1,nX, nY), 'dp_yx is nonzero');
            testCase.verifyEqual(dc_y, zeros(1, nY), 'dc_y is nonzero');
            testCase.verifyEqual(dc_yx, zeros(1,nX, nY), 'dc_yx is nonzero');
            % Check the gradients against analytic results
            testCase.verifyEqual(dp_x, [24.5, -25, -25, -25]*pdf, 'RelTol',1e-12, 'pdf gradient is inaccurate');
            testCase.verifyEqual(dc_x, [-1, 1, 1, 1]*pdf, 'RelTol',1e-12, 'cdf gradient is inaccurate');
            % Check the gradients with respect to the solutions
            pxx_analytic = [550.5, -575, -575, -575;
                            -575, 600, 600, 600;
                            -575, 600, 600, 600;
                            -575, 600, 600, 600] * pdf;
            testCase.verifyEqual(squeeze(dp_xx), pxx_analytic, 'RelTol', 1e-12, 'pdf hessian is inaccurate');
            
            % Verify the first derivatives with respect to x numerically
%             p_fun = @(z) testCase.solver.evalDistribution(z, testCase.P, testCase.w);
%             c_fun = @(z) testCase.outputWrapper2(p_fun, z);
%             dp_est = testCase.finiteDifference(p_fun, testCase.x);
%             dc_est = testCase.finiteDifference(c_fun, testCase.x);
%             testCase.verifyEqual(dp_x, dp_est, 'RelTol',testCase.tol,'pdf gradient does not match numerical gradient to desired precision');
%             testCase.verifyEqual(dc_x ,dc_est, 'RelTol',testCase.tol, 'cdf gradient does not match numerical gradient to desired precision');
%             
            % Verify the second derivatives with respect to x numerically
            dp_fun = @(z) testCase.solver.evalDistribution(z, testCase.P, testCase.w, testCase.dP, testCase.dw);
            px_fun = @(z) testCase.outputWrapper3(dp_fun, z);
            cx_fun = @(z) testCase.outputWrapper4(dp_fun, z);
            pxx_est = testCase.finiteDifference(px_fun, testCase.x);
            cxx_est = testCase.finiteDifference(cx_fun, testCase.x);
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
           % Check the value of the second derivatives
           testCase.verifyEqual(g_xy, [4, 10; 8, 20; -2 6; -2, 6], 'Mixed derivative is inaccurate');
           % Test the value of the hessian using finite differences
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
        %% --- Test for Probabilistic Degeneracy --- $$
        function testDegeneracy(testCase)
           % Import the IsFinite condition
           import matlab.unittest.constraints.IsFinite;
           % Create the degenerate case
           xd = [0,0,1,1]'; % Degenerate because the 1st component is zero
           % Check the standard deviation
           sigma = testCase.solver.ermDeviation(xd, testCase.P, testCase.w);
           testCase.assertEqual(sigma, 0, 'Sigma is nonzero');
           % Check the cost function and it's gradient
           [g, dg] = testCase.solver.ermCost(xd, testCase.P, testCase.w);
           testCase.verifyThat(g, IsFinite, 'cost is not finite or is NaN');
           testCase.verifyThat(dg, IsFinite, 'cost gradient is not finite or is NaN');
           % Check the second derivatives
           [gxx, gxy] = testCase.solver.ermSecondDerivatives(xd, testCase.P, testCase.w, testCase.dP, testCase.dw);
           testCase.verifyThat(gxx, IsFinite, 'cost hessian is not finite or is NaN');
           testCase.verifyThat(gxy, IsFinite, 'cost mixed partials are not finite or are NaN');
           % Compare gradients to numerical finite differneces
%            cost = @(x) testCase.solver.ermCost(x, testCase.P, testCase.w);
%            dg_est = testCase.finiteDifference(cost, xd);
%            testCase.verifyEqual(dg, dg_est, 'RelTol',testCase.tol,'Degenerate gradient does not match numerical gradient to desired precision');
%            dcost = @(x) testCase.solver.costGradients(x, testCase.P, testCase.w, testCase.dP, testCase.dw);
%            ddg_est = testCase.finiteDifference(dcost, xd);
%            testCase.verifyEqual(gxx, ddg_est, 'RelTol',testCase.tol,'Degenerate hessian does not match numerical hessian to desired precision');
        % NOTE: Finite differencing doesn't work in the degenerate case
        % because the backwards difference (x-dx) produces a negative
        % standard deviation. Otherwise, the non-nan values match the
        % derivatives
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
        function f = outputWrapper7(fun, x)
            [~, ~, ~, ~, ~, ~, f] = fun(x);
        end
        function f = outputWrapper8(fun, x)
            [~, ~, ~, ~, ~, ~, ~, f] = fun(x);
        end
    end
end

