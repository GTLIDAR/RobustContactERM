classdef ContactCartOperatorTests < matlab.unittest.TestCase
    %% ContactCartOperatorTests: Unit Tests for Checking the implementation 
    %   of the physical properties of the DifferentiableContactCart.
    %
    %   ContactCartOperatorTests compares the values calculated using
    %   operator notation (DifferentiableContactCart) to those calculated
    %   explicitly (ContactDrivenCart). 
    %
    %   Luke Drnach
    %   December 6, 2019
    
    properties
        explicitModel = ContactDrivenCart();
        operatorModel = DifferentiableContactCart();
        testConfiguration = [2, pi/3, -3*pi/4];
        testVelocity = [-0.1, 0.01, 0.25];
        tol = 1e-12;
    end
    
    methods (Test)
        function compareMassMatrix(testCase)
            % Get the mass matrix and gradient from the explicit model
            [M, dM] = testCase.explicitModel.massMatrix(testCase.testConfiguration);
            % Get the mass matrix and gradient from the operator model
            [M_o, dM_o] = testCase.operatorModel.massMatrix(testCase.testConfiguration);
            
            % Check that they are equal
            testCase.verifyEqual(M, M_o,'RelTol',testCase.tol,'Mass matrices are not equal');
            % Check that the gradients are equal
            testCase.verifyEqual(dM, dM_o,'RelTol',testCase.tol,'Mass gradients are not equal');
        end
        function compareCoriolisMatrix(testCase)
            % Get the coriolis matrix and its gradient
            [C, dC] = testCase.explicitModel.coriolisMatrix(testCase.testConfiguration, testCase.testVelocity);
            % And from the operator model
            [C_o, dC_o] = testCase.operatorModel.coriolisMatrix(testCase.testConfiguration, testCase.testVelocity);
            
            % Check the matrices are equal
            testCase.verifyEqual(C, C_o, 'RelTol',testCase.tol,'Coriolis matrices are not equal');
            % Check the gradients
            testCase.verifyEqual(dC, dC_o, 'RelTol',testCase.tol,'Coriolis gradients are not equal');            
        end
        function compareGravityMatrix(testCase)
            % Get the gravity vectors and their gradients
            [N, dN] = testCase.explicitModel.gravityMatrix(testCase.testConfiguration);
            [N_o, dN_o] = testCase.operatorModel.gravityMatrix(testCase.testConfiguration);
            
            % Check for equality
            testCase.verifyEqual(N, N_o, 'RelTol',testCase.tol,'Gravity vectors are not equal');
            testCase.verifyEqual(dN, dN_o,'RelTol',testCase.tol,'Gravity gradients are not equal');
        end
        function compareJacobianMatrix(testCase)
            % Get the Jacobians and their gradients
            [J, dJ] = testCase.explicitModel.jacobian(testCase.testConfiguration);
            [J_o, dJ_o] = testCase.operatorModel.jacobian(testCase.testConfiguration);
            
            % Check for equality
            testCase.verifyEqual(J, J_o,'RelTol',testCase.tol, 'Jacobians are not equal');
            testCase.verifyEqual(dJ, dJ_o,'RelTol',testCase.tol, 'Jacobian gradients are not equal');
            
        end
    end
    
end

