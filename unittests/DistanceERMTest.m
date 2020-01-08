classdef DistanceERMTest < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
       properties
        solver;
        P;      % LCP Problem Matrix
        w;      % LCP offset vector
        x;      % LCP solution
        dP;     % LCP matrix gradient
        dw;     % LCP offset gradient
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
            %% TESTMEAN: Test the ERMMean function
            
            
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
            %% TESTDEVIATION: Test the ERMDeviation function           
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
                %% --- Tests for the GaussianERM Superclass --- %%
        function testCostFunctionValues(testCase)
            % Get the function value
            [r, dr] = testCase.solver.ermCost(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
            % Check the sizes of the values as a first-pass
            testCase.verifyEqual(size(r), [1,1], 'Residual is the wrong size');
            testCase.verifyEqual(size(dr), [1, numel(testCase.x)],'Gradient is the wrong size');
            % Check for accuracy of the values (?)
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
        end
        function testSecondDerivatives(testCase)
           nX = numel(testCase.x);
           nY = size(testCase.dP, 3);
           % Get the values of the second derivatives of the cost
           [g_xx, g_xy] = testCase.solver.ermSecondDerivatives(testCase.x, testCase.P, testCase.w, testCase.dP, testCase.dw);
           % Check that the sizes are correct
           testCase.verifyEqual(size(g_xx), [nX, nX], 'Hessian g_xx is the wrong size');
           testCase.verifyEqual(size(g_xy), [nX, nY], 'Mixed derivative g_xy is the wrong size');
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
            [f, r] = testCase.solver.solve(testCase.P, testCase.w, testCase.dP, testCase.dw);
            % Check that the solution is the right size
            testCase.verifyEqual(size(f), size(testCase.x));
            testCase.verifyEqual(size(r), [1,1]);
            % Check that the solution is nonnegative
            testCase.verifyTrue(all(f >= 0), 'ERM solution is negative');
        end
    end
    
end