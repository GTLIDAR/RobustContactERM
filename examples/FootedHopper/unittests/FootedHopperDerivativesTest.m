classdef FootedHopperDerivativesTest < matlab.unittest.TestCase
    %FOOTEDHOPPERDERIVATIVESTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % System Properties
        plant = FootedHopper();
        q = [0,3,pi/4, -pi/3, pi/4]';
        dq = [0.1, -0.9, -0.1, 0.5, 0.01]';
        % Finite Differencing Settings
        h = 1e-8;
        tol = 5e-6;
    end
    methods (TestMethodSetup)
        function modifyPlant(testCase)
            testCase.plant.heelFraction = 1/4;
            testCase.plant.lengths(3) = 1/3;
        end
    end
    
    methods (Test)
        function massDerivativesTest(testCase)
          
            % Get the parameters
            q_ = testCase.q;
            [M,dM,HM] = testCase.plant.massMatrix(q_);
            % Check that the mass matrix is the correct size
            testCase.assertEqual(size(M), [numel(q_), numel(q_)], 'Mass matrix is the wrong size');
            % Calculate the derivatives using finite differencing
            dM_est = testCase.finiteDifference(@(x) testCase.plant.massMatrix(x), q_);
            
            hessfun = @(z) testCase.outputWrapper2(@(x) testCase.plant.massMatrix(x),z);
            HM_est = testCase.finiteDifference(@(x) hessfun(x), q_);
            
            % Compare to the implemented derivatives
            testCase.verifyEqual(dM, dM_est, 'RelTol',testCase.tol, 'Mass Matrix gradient is inaccurate');
            testCase.verifyEqual(HM, HM_est, 'RelTol', testCase.tol, 'Mass Matrix Hessian is inaccurate');
            
        end
        function coriolisDerivativesTest(testCase)
            % Get the parameters
            q_ = testCase.q;
            dq_= testCase.dq;
            [C,dC] = testCase.plant.coriolisMatrix(q_,dq_);
            % Check the shape of the coriolis matrix
            testCase.assertEqual(size(C), [numel(q_), numel(q_)], 'Coriolis matrix is the wrong size');
            % Calculate the derivatives using finite differencing
            dC_q = testCase.finiteDifference(@(x) testCase.plant.coriolisMatrix(x, dq_), q_);
            dC_dq = testCase.finiteDifference(@(x) testCase.plant.coriolisMatrix(q_, x), dq_);
            % Compare to implemented derivatives
            testCase.verifyEqual(dC(:,:,1:numel(q_)),dC_q, 'RelTol', testCase.tol, 'Coriolis configuration derivative inaccurate');
            testCase.verifyEqual(dC(:,:,numel(q_)+1:end), dC_dq, 'RelTol', testCase.tol, 'Coriolis velocity derivative inaccurate');
            
        end
        function gravityDerivativesTest(testCase)
            % Get the parameters
            q_ = testCase.q;
            [N,dN] = testCase.plant.gravityMatrix(q_);
            % Check the size
            testCase.assertEqual(size(N), [5,1], 'Gravity vector is the wrong size');
            % Calculate the derivative using finite differencing
            dN_est = testCase.finiteDifference(@(x) testCase.plant.gravityMatrix(x), q_);
            % Compare to the implemented derivatives
            testCase.verifyEqual(dN, dN_est, 'RelTol', testCase.tol, 'Gravity gradient is inaccurate');
            
        end
        function controllerDerivativesTest(testCase)
            % Get the parameters
            q_ = testCase.q;
            [B,dB] = testCase.plant.controllerMatrix(q_);
            % Check teh size
            testCase.assertEqual(size(B),[5,3],'Controller Matrix is the Wrong size');
            % Calculate the derivatives using finite differencing
            dB_est = testCase.finiteDifference(@(x) testCase.plant.controllerMatrix(x), q_);
            
            % Compare to the implemented derivatives
            testCase.verifyEqual(dB, dB_est, 'RelTol', testCase.tol, 'Controller Selector Derivative inaccurate');
            
        end
        function jacobianDerivativesTest(testCase)
            % Get the parameters
            q_ = testCase.q;
            [J, dJ] = testCase.plant.jacobian(q_);
            
            % Check the size
            testCase.assertEqual(size(J), [4,5], 'Jacobian is the wrong size');
            
            % Calculate the derivative using finite differencing
            dJ_est = testCase.finiteDifference(@(x) testCase.plant.jacobian(x), q_);
            
            % Compare to the implemented derivatives
            testCase.verifyEqual(dJ, dJ_est, 'RelTol',testCase.tol, 'Jacobian derivative is inaccurate');
            
        end
        function jacobianTest(testCase)
            q_ = testCase.q;
            J = testCase.plant.jacobian(q_);
            
            J_est = testCase.finiteDifference(@(x)testCase.plant.kinematics(x), q_);
            J_est = reshape(J_est, [4,5]);
            testCase.verifyEqual(J, J_est,'RelTol',testCase.tol,'Jacobian is inaccurate');
            
        end
    end
    methods
        
        function df = finiteDifference(obj, fun, x)
            
            f = fun(x);
            
            df =  cell(1, numel(x));
            % Shifts
            dx = obj.h/2 * eye(numel(x));
            % Central differencing
            for n = 1:numel(x)
                f1 = fun(x + dx(:,n));
                f2 = fun(x - dx(:,n));
                df{n} = (f1 - f2)./obj.h;
            end
            % Collect the outputs together
            s = size(f);
            if s(2) == 1
                dim = 2;
            else
                dim = numel(s) + 1;
            end
            df = cat(dim, df{:});
        end
    end
    methods (Static)
        function fun2 = outputWrapper2(fun, x)
            [~, fun2] = fun(x);
        end
    end
end

