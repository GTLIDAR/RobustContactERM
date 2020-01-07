classdef FrictionERMTest < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Luke Drnach
    %   January 7, 2019
    
    properties
        solver;
        A;      % LCP Problem Matrix
        b;      % LCP offset vector
        z;      % LCP solution
        dA;     % LCP matrix gradient
        db;     % LCP offset gradient
    end
    
    methods (TestMethodSetup)
        function createExample(testCase)
            % Create the plant needed to create a solver
            plant = Block();    % Use a block with friction as our example
            mu = 0.5;           % Test value for friction coefficient
            sig = 0.1;          % Test value for friction uncertainty
            testCase.solver = UncertainFrictionERM(plant,mu,sig);
            % Example LCP parameters
            testCase.A = [1    2  0 0;      %LCP problem matrix
                          2    1  0 0;
                          0    0  1 1;
                          0.5 -1 -1 0];
            testCase.b = [-2, 2, -2, 0]';   %LCP problem vector
            testCase.z = [2, 0, 1, 1]';     %LCP Solution vector
            % Example LCP Gradients
            testCase.dA = zeros(4,4,2);
            testCase.dA(:,:,1) = [1 -1 0 0; -1 2 1 0; 0 1 0 0; 0 0 0 0];
            testCase.dA(:,:,2) = [1 0 1 0; 0 0 -1 0; 1 -1 2 0; 0 0 0 0];
            testCase.db = [0 2; 1 1; -1 -1; 0 0]; 
        end
    end
    methods (Test)
        function testMean(testCase)
            %% TESTMEAN: Test the ERMMean function
            
            % Get the problem parameters
            P = testCase.A;
            w = testCase.b;
            x = testCase.z;
            dP = testCase.dA;
            dw = testCase.db;
            
            % Get the values for the mean
            [mu, dmu_x, dmu_y, dmu_xx, dmu_xy] = testCase.solver.ermMean(x,P,w,dP,dw);
            
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
            
            % Get the problem parameters
            P = testCase.A;
            w = testCase.b;
            x = testCase.z;
            dP = testCase.dA;
            dw = testCase.db;
            
            % Get the value of the standard deviation
            [m_sigma, dsigma_x, dsigma_y, dsigma_xx, dsigma_xy]  = testCase.solver.ermDeviation(x, P, w, dP, dw);
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
    end
end

