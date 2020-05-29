classdef SemiImplicitTest < matlab.unittest.TestCase
    % Numerical tests of the values and gradients of the Semi-Implicit
    % Integration Method in SemiImplicitTrajectoryOptimization
    
    
    
    properties
        block;
        h = 1e-8;
        tol = 1e-6;
    end
    
    methods (TestMethodSetup)
        function createExamples(testCase)
            % Create the pushblock example
            plant = Block();
            plant.timestep = 0.01;
            x0 = [0, plant.height/2, 0, 0]';
            x1 = [1, plant.height/2, 0, 0]';
            u = 1;
            Tf = 1;
            N = 101;
            prob = SemiImplicitTrajectoryOptimization(plant, N, Tf);
            testCase.block = struct('plant',plant,'x0',x0,'x1',x1,'u0',u,'Tf',Tf,'N',N,...
                'solver',prob);
            
        end
    end
    
    methods (Test)
        function pushBlockTest_NoContact(testCase)
            example = testCase.block;
            x0 = example.x0;
            x1 = example.x1;
            x0(2) = 2;
            x1(2) = 2;
            [~, df] = example.solver.semiimplicit_constraint_fun(0, x0, x1, example.u0);
            
            % Numerical finite difference check for h
            fun = @(z) example.solver.semiimplicit_constraint_fun(z, x0, x1, example.u0);
            df_est = testCase.finiteDifference(fun, 0);
            testCase.verifyEqual(df(:,1),df_est,'RelTol',testCase.tol,'Time Derivative Inaccurate');
            % Numerical finite difference check for x0
            fun = @(z) example.solver.semiimplicit_constraint_fun(0, z, x1, example.u0);
            df_est = testCase.finiteDifference(fun, example.x0);
            testCase.verifyEqual(df(:,2:5),df_est, 'RelTol',testCase.tol, 'State Derivative at 0 Inaccurate');
            % Numerical finite difference check for x1
            fun = @(z) example.solver.semiimplicit_constraint_fun(0, x0, z, example.u0);
            df_est = testCase.finiteDifference(fun, example.x1);
            testCase.verifyEqual(df(:,6:9), df_est, 'RelTol',testCase.tol,'State Derivative at 1 Inaccurate');
            % Numerical finite difference check for u0
            fun = @(z) example.solver.semiimplicit_constraint_fun(0, x0, x1, z);
            df_est = testCase.finiteDifference(fun, example.u0);
            testCase.verifyEqual(df(:,10),df_est, 'RelTol',testCase.tol,'Control Derivative Inaccurate');
        end
        function pushBlockTest_Contact(testCase)
           example = testCase.block;
           [f, df] = example.solver.semiimplicit_constraint_fun(0, example.x0, example.x1, example.u0);
           
           % Numerical finite difference check for h
           fun = @(z) example.solver.semiimplicit_constraint_fun(z, example.x0, example.x1, example.u0);
           df_est = testCase.finiteDifference(fun, 0);
           testCase.verifyEqual(df(:,1),df_est,'RelTol',testCase.tol,'Time Derivative Inaccurate');
           % Numerical finite difference check for x0
           fun = @(z) example.solver.semiimplicit_constraint_fun(0, z, example.x1, example.u0);
           df_est = testCase.finiteDifference(fun, example.x0);
           testCase.verifyEqual(df(:,2:5),df_est, 'RelTol',testCase.tol, 'State Derivative at 0 Inaccurate');
           % Numerical finite difference check for x1
           fun = @(z) example.solver.semiimplicit_constraint_fun(0, example.x0, z, example.u0);
           df_est = testCase.finiteDifference(fun, example.x1);
           testCase.verifyEqual(df(:,6:9), df_est, 'RelTol',testCase.tol,'State Derivative at 1 Inaccurate');
           % Numerical finite difference check for u0
           fun = @(z) example.solver.semiimplicit_constraint_fun(0, example.x0, example.x1, z);
           df_est = testCase.finiteDifference(fun, example.u0);
           testCase.verifyEqual(df(:,10),df_est, 'RelTol',testCase.tol,'Control Derivative Inaccurate');
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

