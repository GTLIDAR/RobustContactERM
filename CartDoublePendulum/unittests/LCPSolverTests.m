classdef LCPSolverTests < matlab.unittest.TestCase
   
    properties
       solver = PathLCPSolver();
       example = [];
       tolerance = 1e-8;
    end
    
    methods (TestClassSetup)
        function createExample(testCase)
            % The example LCP problem
            testCase.example.A = [2, 1; 1, 2];
            testCase.example.b = [-2, 1]';
            testCase.example.x = [1; 0];
            % The example derivative problem
            testCase.example.dA = zeros(2,2,2);
            testCase.example.dA(:,:,1) = [1 ,-1; 0, 1];
            testCase.example.dA(:,:,2) = [-1, 0; 2, 1];
            testCase.example.db = [1, -1; 0, 2];
            %testCase.example.dx = -[4, -8; -2, 10]./3; % This is the gradient if we don't make the gradient 0 when x is 0
            testCase.example.dx = -[4, -8; 0, 0]./3;
        end
    end
    
    methods (Test)
        function testSolver(testCase)
           % Solve the LCP problem
           y = testCase.solver.solve(testCase.example.A, testCase.example.b);
           % Compare the solution to the example solution
           testCase.verifyEqual(y, testCase.example.x, 'RelTol',testCase.tolerance);
        end
        function testGradient(testCase)
            % Calculate the gradient for the problem
            dx = testCase.solver.gradient(testCase.example.x, testCase.example.A,...
                testCase.example.b, testCase.example.dA, testCase.example.db);
            % Compare to the analytical gradient
            testCase.verifyEqual(dx,testCase.example.dx,'RelTol',testCase.tolerance);
        end
        function testSolverAndGradient(testCase)
            % Solve the LCP problem
            y = testCase.solver.solve(testCase.example.A, testCase.example.b);
            % calculate the gradient using the LCP solution
            dx = testCase.solver.gradient(y, testCase.example.A,...
                testCase.example.b, testCase.example.dA, testCase.example.db);
            % Compare to the analytical gradient
            testCase.verifyEqual(dx,testCase.example.dx,'RelTol',testCase.tolerance);
        end
    end
    
end