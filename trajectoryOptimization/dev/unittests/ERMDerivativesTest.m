classdef ERMDerivativesTest < matlab.unittest.TestCase
    %ERMDERIVATIVESTEST Summary of this class goes here

    %   Luke Drnach
    %   March 12, 2020
    
    properties
        optimizer;
        inContact;
        noContact;
        ermTestVals;
        % Finite Difference parameters
        h = 1e-6;
        tol = 1e-6;
    end
    
    methods (TestMethodSetup)
        function createExamples(testCase)
           % Create teh plant model
           plant = ContactDrivenCart;
           plant.timestep = 0.01;
           plant.cartHeight = 1.5;
           
           % Specify in-contact and no-contact conditions
           x0 = [0,0]'; % Endpoint in contact;
           xf = [0,1]'; % Endpoint not in contact;
           % Calculate corresponding states
           q0 = plant.inverseKinematics(x0);
           qf = plant.inverseKinematics(xf);
           % Set velocities equal to 0
           x0 = [q0; zeros(3,1)];
           xf = [qf; zeros(3,1)];
           
           % Store values
           testCase.inContact = struct('x',x0,'h',plant.timestep,'f',200);
           testCase.noContact = struct('x',xf,'h',plant.timestep,'f',0);
           
           % Create optimization problem
           testCase.optimizer = RobustContactImplicitTrajectoryOptimizer(plant, 101, 1);
           
           % Create some other test values
           testCase.ermTestVals = struct('x',[1,3,0,0.5,0]','mu',[0,0.5,1,0.1,0.2]', 'sigma',[1,0.4,1.2,0,0]');
        end
    end
    methods (Test)
        function testERMValues(testCase)
            
           % Pull the test values from the testCase
           x_ = testCase.ermTestVals.x;
           mu_ = testCase.ermTestVals.mu;
           sigma_ = testCase.ermTestVals.sigma;
           
           % Get the value of the ERM from the old implementation
           [f, df_x, df_mu, df_sigma] = testCase.ermCostGaussian_old(x_, mu_, sigma_);
           df = [df_x, df_mu, df_sigma];
           % Get the values of the ERM from the current implementation
           [g, dg] = testCase.optimizer.ermCostGaussian(x_, mu_, sigma_);
           
           % Compare the two values
           testCase.verifyEqual(f, g, 'AbsTol',1e-12, 'Current implementation does not match old implementation');
           % Compare the two gradients
           testCase.verifyEqual(df(:,1:numel(x_)), dg(:,1:numel(x_)), 'AbsTol',1e-12, 'Current gradient implementation does not match old implementation in X');
           testCase.verifyEqual(df(:,numel(x_)+1:2*numel(x_)), dg(:,numel(x_)+1:2*numel(x_)), 'AbsTol',1e-12, 'Current gradient implementation does not match old implementation in MU');
           testCase.verifyEqual(df(:,2*numel(x_)+1:end), dg(:,2*numel(x_)+1:end), 'AbsTol',1e-12, 'Current gradient implementation does not match old implementation in SIGMA');
        end
        function testERMGradient(testCase)
            
            % Pull the test values from the testCase
            x_ = testCase.ermTestVals.x;
            mu_ = testCase.ermTestVals.mu;
            sigma_ = testCase.ermTestVals.sigma;
            nX = numel(x_);
            
            % Calculate the ERM Gradient
            [~, df] = testCase.optimizer.ermCostGaussian(x_, mu_, sigma_);
            
            % Calculate the Gradient using finite differencing
            df_x = testCase.finiteDifference(@(z) testCase.optimizer.ermCostGaussian(z, mu_, sigma_), x_);
            testCase.verifyEqual(df(:,1:nX), df_x, 'AbsTol',testCase.tol, 'Derivative wrt X inaccurate');
            
            % Gradient wrt the mean
            df_mu = testCase.finiteDifference(@(z) testCase.optimizer.ermCostGaussian(x_, z, sigma_), mu_);
            testCase.verifyEqual(df(:,nX+1:2*nX), df_mu, 'AbsTol', testCase.tol, 'Derivative wrt MU inaccurate');
            
            % Gradient wrt the variance
            df_sigma = testCase.forwardDifference(@(z) testCase.optimizer.ermCostGaussian(x_, mu_, z), sigma_);
            testCase.verifyEqual(df(:,2*nX+1:end), df_sigma, 'AbsTol', testCase.tol, 'Derivative wrt SIGMA inaccurate');
        end
        function testERMHessian(testCase)
            
            % Pull the test values from the testCase
            x_ = testCase.ermTestVals.x;
            mu_ = testCase.ermTestVals.mu;
            sigma_ = testCase.ermTestVals.sigma;
            nX = numel(x_);
            
            % Calculate the ERM Hessian
            [~, ~, Hf] = testCase.optimizer.ermCostGaussian(x_, mu_, sigma_);
            
            hessfun = @(x, y, z)testCase.outputWrapper2(@(x,y,z)testCase.optimizer.ermCostGaussian(x,y,z), x, y, z);
            
            % Check the hessian wrt x
            Hf_x = testCase.finiteDifference(@(z)hessfun(z, mu_, sigma_), x_);
            testCase.verifyEqual(Hf(:,:,1:nX), Hf_x, 'AbsTol',testCase.tol, 'Hessian wrt X inaccurate');
            
            % Check the hessian wrt mu
            Hf_mu = testCase.finiteDifference(@(z)hessfun(x_, z, sigma_), mu_);
            testCase.verifyEqual(Hf(:,:,nX+1:2*nX), Hf_mu, 'AbsTol',testCase.tol, 'Hessian wrt MU inaccurate');
            
            % Check the hessian wrt sigma
            Hf_sigma = testCase.forwardDifference(@(z)hessfun(x_, mu_, z), sigma_);
            testCase.verifyEqual(Hf(:,:,2*nX+1:end), Hf_sigma, 'AbsTol', testCase.tol, 'Hessian wrt SIGMA inaccurate');            
        end
        function testDistanceGradientContact(testCase)
            
            h_ = testCase.inContact.h;
            x_ = testCase.inContact.x;
            l_ = testCase.inContact.f;
            [~, df] = testCase.optimizer.normalDistanceERMCost(h_, x_, l_);
            % Test the gradient wrt timestep
            df_h = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(z, x_, l_), h_);
            testCase.verifyEqual(df(:,1), df_h, 'AbsTol',testCase.tol, 'Derivative wrt time step inaccurate');
            % Test the gradient wrt state
            df_x = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, z, l_), x_);
            testCase.verifyEqual(df(:,2:1+numel(x_)), df_x, 'AbsTol',testCase.tol, 'Derivative wrt state inaccurate');
            % Test the gradient wrt force
            df_l = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, x_, z), l_);
            testCase.verifyEqual(df(:,2+numel(x_):end), df_l, 'AbsTol', testCase.tol, 'Derivative wrt force inaccurate');
        end
        function testDistanceGradientNoContact(testCase)
            
            h_ = testCase.noContact.h;
            x_ = testCase.noContact.x;
            l_ = testCase.noContact.f;
            [~, df] = testCase.optimizer.normalDistanceERMCost(h_, x_, l_);
            % Test the gradient wrt timestep
            df_h = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(z, x_, l_), h_);
            testCase.verifyEqual(df(:,1), df_h, 'AbsTol',testCase.tol, 'Derivative wrt time step inaccurate');
            % Test the gradient wrt state
            df_x = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, z, l_), x_);
            testCase.verifyEqual(df(:,2:1+numel(x_)), df_x, 'AbsTol',testCase.tol, 'Derivative wrt state inaccurate');
            % Test the gradient wrt force
            df_l = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMCost(h_, x_, z), l_);
            testCase.verifyEqual(df(:,2+numel(x_):end), df_l, 'AbsTol', testCase.tol, 'Derivative wrt force inaccurate');
        end
        function testDistanceConstraintHessianContact(testCase)
           
            h_ = testCase.inContact.h;
            x_ = testCase.inContact.x;
            l_ = testCase.inContact.f;
            [~, df] = testCase.optimizer.normalDistanceERMGradientConstraint([h_; x_; l_]);
            % Test the gradient wrt timestep
            df_h = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMGradientConstraint([z; x_; l_]), h_);
            testCase.verifyEqual(df(:,1), df_h, 'AbsTol',testCase.tol, 'Derivative wrt time step inaccurate');
            % Test the gradient wrt state
            df_x = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMGradientConstraint([h_; z; l_]), x_);
            testCase.verifyEqual(df(:,2:1+numel(x_)), df_x, 'AbsTol',testCase.tol, 'Derivative wrt state inaccurate');
            % Test the gradient wrt force
            df_l = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMGradientConstraint([h_; x_; z]), l_);
            testCase.verifyEqual(df(:,2+numel(x_):end), df_l, 'AbsTol', testCase.tol, 'Derivative wrt force inaccurate');
        end
        function testDistanceConstraintHessianNoContact(testCase)
            
            
            h_ = testCase.noContact.h;
            x_ = testCase.noContact.x;
            l_ = testCase.noContact.f;
            [~, df] = testCase.optimizer.normalDistanceERMGradientConstraint([h_; x_; l_]);
            % Test the gradient wrt timestep
            df_h = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMGradientConstraint([z; x_; l_]), h_);
            testCase.verifyEqual(df(:,1), df_h, 'AbsTol',testCase.tol, 'Derivative wrt time step inaccurate');
            % Test the gradient wrt state
            df_x = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMGradientConstraint([h_; z; l_]), x_);
            testCase.verifyEqual(df(:,2:1+numel(x_)), df_x, 'AbsTol',testCase.tol, 'Derivative wrt state inaccurate');
            % Test the gradient wrt force
            df_l = testCase.finiteDifference(@(z) testCase.optimizer.normalDistanceERMGradientConstraint([h_; x_; z]), l_);
            testCase.verifyEqual(df(:,2+numel(x_):end), df_l, 'AbsTol', testCase.tol, 'Derivative wrt force inaccurate');
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
        function df = forwardDifference(obj, fun, x)
            
            f = fun(x);
            
            df =  cell(1, numel(x));
            % Shifts
            dx = obj.h * eye(numel(x));
            % Central differencing
            for n = 1:numel(x)
                f1 = fun(x + dx(:,n));
                df{n} = (f1 - f)./obj.h;
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
        function fun2 = outputWrapper2(fun, varargin)
            [~, fun2] = fun(varargin{:});
        end
        function [f, dfx, dfmu, dfsigma, d2fx, d2fmu, d2fsigma] = ermCostGaussian_old(x, mu, sigma)
            %% ERMCOSTGassian: Cost function for the Expected Residual Minimization with Gaussian variables
            %
            %
            
            % Initialize the outputs
            f = zeros(length(x),1);
            dfx = zeros(length(x),1);
            dfmu = dfx;
            dfsigma = dfx;
            d2fx  = zeros(length(x), 1);
            d2fmu = d2fx;
            d2fsigma = d2fx;
            % Filter out any degenerate distributions (no variance cases)
            degenerate = (sigma == 0);
            % Save the degenerate means for the limiting case
            mu_degenerate = mu(degenerate);
            x_degenerate = x(degenerate);
            % Filter out the remaining variables
            x = x(~degenerate);
            sigma = sigma(~degenerate);
            mu = mu(~degenerate);
            if any(~degenerate)
                % Calculate the pdf and cdf values
                pdf = normpdf(x, mu, sigma);
                cdf = normcdf(x, mu, sigma);
                % The ERM cost function for Gaussian variables with nonzero
                % variance
                f(~degenerate) = x.^2 - sigma.^2 .* (x + mu) .* pdf + (sigma.^2 + mu.^2 - x.^2) .* cdf;
                % Calculate the derivaties of pdf and cdf wrt x, mu, and sigma
                tau = (x - mu)./sigma;
                dfx(~degenerate) = 2*x.*(1 - cdf);
                dfmu(~degenerate) = 2.*(mu.*cdf - sigma.^2 .* pdf);
                dfsigma(~degenerate) = 2.*sigma .*(cdf - x.*pdf);
                % The Hessian Values
                d2fx(~degenerate) = 2.*(1 - x.*pdf - cdf);
                d2fmu(~degenerate) = 2.*(cdf - x.*pdf);
                d2fsigma(~degenerate) = -2.*((x - mu) + tau.^2).*pdf + 2.*cdf;
            end
            if any(degenerate)
                % Handle the degenerate case
                x_idx = (x_degenerate < mu_degenerate);
                % Initialize holding arrays
                g = zeros(length(x_degenerate), 1);
                dg_x = g;
                dg_mu = g;
                d2g_x = g;
                d2g_mu = g;
                d2g_sigma = g;
                % Calculate function values for the degenerate case
                g(x_idx) = x_degenerate(x_idx).^2;
                g(~x_idx) = mu_degenerate(~x_idx).^2;
                % Calculate gradients
                dg_x(x_idx) = 2.*x_degenerate(x_idx);
                dg_mu(~x_idx) = 2.*mu_degenerate(~x_idx);
                % Hessians
                d2g_x(x_idx) = 2;
                d2g_mu(~x_idx) = 2;
                d2g_sigma(~x_idx) = 2;
                % Map back to the output
                f(degenerate) = g;
                dfx(degenerate) = dg_x;
                dfmu(degenerate) = dg_mu;
                dfsigma(degenerate) = 0;
                d2fx(degenerate) = d2g_x;
                d2fmu(degenerate) = d2g_mu;
                d2fsigma(degenerate) = d2g_sigma;
            end
            % Convert the gradients into vectors/matrices
            dfx = diag(dfx);
            dfmu = diag(dfmu);
            dfsigma = diag(dfsigma);
            % Convert the hessians into matrices/3D arrays
            I = eye(numel(f)) .* reshape(eye(numel(f)), numel(f), 1, numel(f)); % Useful for dealing out numbers to the right positions
            d2fx = I.*d2fx;
            d2fmu = I.*d2fmu;
            d2fsigma = I.*d2fsigma;
        end
    end
end

