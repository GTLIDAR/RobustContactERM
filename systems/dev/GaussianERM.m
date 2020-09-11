classdef GaussianERM < ContactSolver
    %% GAUSSIANERM: Implements template for the Expected Residual Minimization approach to solving SLCP with Normally distributed uncertianty
    %
    %   GaussianERM is an abstract class that implements common
    %   functionality to all versions of the Expected Residual Minimization
    %   (ERM) technique to solving Stochastic Linear Complementarity
    %   Problems (SLCP) with Gaussian distributed problem variables.
    %   GaussianERM solves the SLCP problem using FMINCON, and can also
    %   return the gradient of the solution with respect to other,
    %   domain-specific variables. 
    %
    %   To implement a GaussianERM, the user must provide functions for
    %   determining the mean and variance of the Gaussian Distribution of
    %   the LCP slack variable z. 
    %
    %   See also: UncertainFrictionERM, UncertainDistanceERM
    %
    %   Luke Drnach
    %   December 18, 2019.
    
    % TODO: Check gradient implementation 
    %       Warm start ERM Solver - Eliminate problems with singularities
    
    properties
        mu;              % Mean of the distribution over the uncertain parameter (NOT the mean used by the final ERM formula)
        sigma;           % Standard deviation of the distribution over the uncertain parameters
        options;         % Options structure for FMINUNC
        guess = [];      % Initial guess for ERM
        Reg = [];        % Regularizer matrix
    end
    properties (SetAccess = protected)
        uncertainIdx;    %Logical index indicating which variables are affected by the uncertainty
    end
   methods (Abstract)
       [m_mu, dmu_x, dmu_xx, dmu_y, dmu_xy] = ermMean(obj, x, P, w, varargin);
       [m_sigma, dsigma_x, dsigma_xx, dsigma_y, dsigma_xy] = ermDeviation(obj, x, P, w, varargin);
   end
   
   methods
       function obj = GaussianERM(m, s, idx)
           %% GaussianERM: Creates an instance of the GaussianERM class
           %
           %   Arguments:
           %       m:      Ns x 1 double, the constant means of the
           %               Gaussian distribution describing the parameter
           %               uncertainty.
           %       s:      Ns x 1 double, the standard deviations of the
           %               Gaussian distribution describing the parameter
           %               uncertainty.
           %       idx:    N x 1 logical, where true indicates the
           %               variables subject to uncertainty.
           %   Return Values:
           %       OBJ:     An instance of the GaussianERM class
           %
           %   Notes: Ns is the number of variables subject to uncertainty.
           %   In general, the number of TRUE values in IDX and Ns should
           %   be equal.
           
           obj.mu = m;
           obj.sigma = s;
           obj.uncertainIdx = idx;
           obj.options = optimoptions('fmincon','Algorithm','interior-point',...
                'SpecifyObjectiveGradient',true,'display','none');
           obj.Reg = zeros(numel(idx));
       end
       function [f, r] = solve(obj, P, w)
           %% SOLVE: Calculate the solution to the ERM problem
           %
           %   SOLVE uses unconstrained optimization to calculate the
           %   solution to the ERM problem.
           %
           %   Arguments:
           %       OBJ:    an instance of the GaussianERM class
           %       P:      NxN double, the LCP problem matrix
           %       w:      Nx1 double, the LCP problem vector
           %
           %   Return Values:
           %       f:      Nx1 double, the solution the SLCP-ERM problem
           %       r:      scalar double, the residual (cost function
           %               evaluation)
           
           % Calculate an initial guess
           %f0 = pathlcp(P, w);
           f0 = ones(numel(w), 1);
           m_mu = obj.ermMean(f0, P, w);
           f0(obj.uncertainIdx,:) = m_mu;   % Set stochastic variables to the mean value
           % Cost function
           cost = @(x) ermCost(obj, x, P, w);
           % Hessian function
           obj.options = optimoptions(obj.options, 'HessianFcn',@(x, lambda) ermSecondDerivatives(obj, x, P, w));
           % Solve the problem
           [f, r, exitflag] = fmincon(cost, f0, -P, w, [], [], zeros(size(f0)), [], [], obj.options);
           if exitflag<=0
              warning('FMINCON might have terminated before a solution was found'); 
           end
           obj.guess =f0;
       end

       function df = gradient(obj, x, P, w, dP, dw)
           %% ERM_COST: Objective function for the ERM approach
           %
           %   ERM_COST returns the sum of the expected residuals for the
           %   Stochastic LCP problem. Minimizing the sum of residuals
           %   produces a solution to the SLCP problem.
           %
           %   Arguments:
           %       OBJ:    an instance of the GaussianERM class
           %       x:      Nx1 double, the current guess of the LCP solution
           %       P:      NxN double, the LCP problem matrix
           %       w:      Nx1 double, the LCP problem vector
           %       dP:     NxNxM double, array of derivatives of P
           %       dw:     NxM double, array of derivatives of w
           %
           %   * For contact problems, M is usually the number of state
           %   variables + the number of control variables.
           %
           %   Return Values:
           %       df:    NxM double, the derivatives of the solution f
           %              with respect to the other problem variables (usually state and control)
           % Get the second derivatives
           [g_xx, g_xy] = obj.ermSecondDerivatives(x, P, w, dP, dw);
           % Solve the system of equations to get the gradient
           df = - g_xx\g_xy;
           % NOTE: Check if we need to zero out the components of the
           % gradient corresponding to the zero solution.
       end
       function [r, dr] = ermCost(obj, x, P, w)
           %% ERMCOST: Objective function for the ERM approach
           %
           %   ERMCOST returns the sum of the expected residuals for the
           %   Stochastic LCP problem. Minimizing the sum of residuals
           %   produces a solution to the SLCP problem.
           %
           %   Arguments:
           %       OBJ:    an instance of the GaussianERM class
           %       x:      Nx1 double, the current guess of the LCP solution
           %       P:      NxN double, the LCP problem matrix
           %       w:      Nx1 double, the LCP problem vector
           %
           %   * For contact problems, M is usually the number of state
           %   variables + the number of control variables.
           %
           %   Return Values:
           %       r:     scalar double, the sum of expected residuals
           %       dr:    Nx1 double, partial derivative of the residual
           %              with respect to the LCP parameter x
           
           % Calculate the NCP cost for the entire problem
           z = P * x + w;
           r = min(x, z).^2;
           % Get the mean and standard deviation for each problem
           [m_mu, dmu_x,] = obj.ermMean(x, P, w);
           [m_sigma, dsigma_x] = obj.ermDeviation(x, P, w);
           % Re-calculate the cost for the stochastic variables
           x_s = x(obj.uncertainIdx,:);
           % Calculate the ERM residual
           [pdf, cdf, dp_x, dc_x] = obj.evalDistribution(x, P, w);
           % Check for the degenerate case
           %zeroIdx = (m_sigma == 0);
           % Calculate the residuals for only the stochastic part
           r_s = x_s.^2 - m_sigma.^2 .* (x_s + m_mu).*pdf + ...
               (m_sigma.^2 + m_mu.^2 - x_s^2).*cdf;
           % The degenerate distribution case
           %r_s(zeroIdx,:) = m_mu(zeroIdx).^2;
           % Combine stochastic and deterministic residuals
           r(obj.uncertainIdx,:) = r_s;
           % Sum all the residuals together
           r = sum(r);
           % Add the regularizer to the cost
           r = r + (x' * obj.Reg * x)/2;

           %% Calculate the gradient
           if nargout == 2
               nX = numel(x);
               % First derivative for deterministic variables
               dr = 2 * z .* P;
               Dx = diag(x);
               dr(x < z,:) = 2.*Dx(x < z,:);
               % First derivative for stochastic variables
               % Get the stochastic variables
               % Get the uncertain variables
               x_s = x(obj.uncertainIdx,:);
               del_nj = zeros(length(x_s),nX);
               del_nj(:,obj.uncertainIdx) = eye(sum(obj.uncertainIdx));
               % Calculate the terms multiplying the distribution values
               pdf_common = 2.*m_sigma.*dsigma_x .* (x_s + m_mu) + m_sigma.^2 .*(del_nj + dmu_x);
               dp_x_common = m_sigma.^2 .*(x_s + m_mu);
               cdf_common = 2*(m_sigma.*dsigma_x + m_mu.*dmu_x - x_s.*del_nj);
               dc_x_common = m_sigma.^2 + m_mu.^2 - x_s.^2;
               % The derivative of the stochastic cost
               df_x = 2.*Dx(obj.uncertainIdx,:) - pdf_common .* pdf - dp_x_common .* dp_x + ...
                   cdf_common .* cdf + dc_x_common .* dc_x;
               % Degenerate case
               %df_x(zeroIdx,:) = 2.*m_mu(zeroIdx).*dmu_x(zeroIdx,:);
               % Combine the two costs
               dr(obj.uncertainIdx,:) = df_x;
               dr = sum(dr, 1);
               % Add the regularizer to the gradient
               dr = dr + x'*obj.Reg;
           end
       end
   end
   methods (Sealed)
       function [dg_xx, dg_xy] = ermSecondDerivatives(obj, x, P, w, dP, dw)
           %% ERMSECONDDERIVATIVES: Hessian and Mixed Partial Derivatives of the Cost
           %
           %   ermSecondDerivatives returns Hessian dg_xx and the mixed
           %   partial derivatives dg_yx such that the solution gradient dx_y
           %   satisfies:
           %            dg_xx * dx_y = -dg_yx
           %
           %   Arguments:
           %       OBJ:    an instance of the GaussianERM class
           %       x:      Nx1 double, decision variables of the ERM
           %               problem
           %       P:      NxN double, the LCP problem matrix
           %       w:      Nx1 double, the LCP problem vector
           %       dP:     NxNxM double, array of derivatives of P
           %       dw:     NxM double, array of derivatives of w
           %
           %   * For contact problems, M is usually the number of state
           %   variables + the number of control variables.
           %
           %   Return Values:
           %       g_xx:  NxN double, a matrix of second derivatives of the
           %              cost function g with respect to the decision
           %              variable x
           %       g_yx:  NxM double, a matrix of mixed partial derivatives
           %              of the cost function g with respect to the 
           %              decision variables x and the parameters y
           
           %% Setup
           % Calculate the distribution values
           if nargout == 2
               [pdf, cdf, dp_x, dc_x, dp_xx, dc_xx, dp_y, dc_y, dp_xy, dc_xy] = obj.evalDistribution(x, P, w, dP, dw);
               [m_sigma, dsigma_x, dsigma_xx, dsigma_y, dsigma_xy]  = obj.ermDeviation(x, P, w, dP, dw);
               [m_mu, dmu_x, dmu_xx, dmu_y, dmu_xy] = obj.ermMean(x, P, w, dP, dw);
           else
               % Calculate the distribution values
               [pdf, cdf, dp_x, dc_x, dp_xx, dc_xx] = obj.evalDistribution(x, P, w);
               [m_sigma, dsigma_x, dsigma_xx]  = obj.ermDeviation(x, P, w);
               [m_mu, dmu_x, dmu_xx] = obj.ermMean(x, P, w);
           end
           % Check for degenerate distributions
           %zeroIdx = (m_sigma == 0);
           % Reshape some of the variables to make the multiplications work
           dsigma_xk = reshape(dsigma_x,size(dsigma_x,1),1,size(dsigma_x,2));
           dmu_xk = reshape(dmu_x,size(dmu_x,1),1,size(dmu_x,2));
           dp_xk = reshape(dp_x,size(dp_x,1),1,size(dp_x,2));
           dc_xk = reshape(dc_x,size(dc_x,1),1,size(dc_x,2));
           
           % Calculate the slack variables
           z = P*x + w;
           
           nX = numel(x);
           %% Hessian
           % Terms where z < x
           P2 = reshape(P,[nX,1,nX]);   % Broadcast P to make multiplication easier
           dg_xx = 2 * P .* P2;
           % Terms where x < z
           I = eye(nX);
           delta_3 = I .* reshape(I,[nX,1,nX]);
           dg_xx(x < z,:,:) = 2*delta_3(x < z,:,:);
           % Uncertain terms
           % Get the uncertain variables
           x_s = x(obj.uncertainIdx,:);
           % Delta functions
           del_nj = zeros(length(x_s),nX);
           del_nj(:,obj.uncertainIdx) = eye(sum(obj.uncertainIdx));
           del_nk = reshape(del_nj,size(del_nj,1),1,size(del_nj,2));
           
           % Calculate terms common to both the hessian and the mixed
           % derivative
           pdf_common = 2.*m_sigma.*dsigma_x .* (x_s + m_mu) + m_sigma.^2 .*(del_nj + dmu_x);
           dp_x_common = m_sigma.^2 .*(x_s + m_mu);
           cdf_common = 2*(m_sigma.*dsigma_x + m_mu.*dmu_x - x_s.*del_nj);
           dc_x_common = m_sigma.^2 + m_mu.^2 - x_s.^2;
           % Calculate the terms that multiply each of the distribution
           % values for the hessian
           pdf_mult_x = 2*(dsigma_x.*dsigma_xk + m_sigma.*dsigma_xx).*(x_s + m_mu) + 2 .* m_sigma.*dsigma_x.*(del_nk + dmu_xk) + ...
               2.*m_sigma.*dsigma_xk .*(del_nj + dmu_x) + m_sigma.^2.*dmu_xx;
           dp_x_mult_x = 2.*m_sigma.*dsigma_xk.*(x_s+m_mu) + m_sigma.^2.*(del_nk + dmu_xk);
           cdf_mult_x = 2.*(dsigma_xk.*dsigma_x + m_sigma.*dsigma_xx + dmu_xk.*dmu_x + m_mu.*dmu_xx - del_nk.*del_nj);
           dc_x_mult_x = 2.*(m_sigma.*dsigma_xk + m_mu.*dmu_xk - x_s.*del_nk);
           % Calculate the hessian
           df_xx = 2.*del_nk.*del_nj - pdf_mult_x .* pdf - pdf_common .* dp_xk - ...
               dp_x_mult_x .* dp_x - dp_x_common .* dp_xx + cdf_mult_x .* cdf + cdf_common .* dc_xk + ...
               dc_x_mult_x .* dc_x + dc_x_common .* dc_xx;
           % Correct for degeneracy
           %df_xx(zeroIdx,:,:) = 2.*(dmu_xk(zeroIdx,:,:).*dmu_x(zeroIdx,:) + m_mu(zeroIdx).*dmu_xx(zeroIdx,:,:));
           % Combine the partial derivatives together
           dg_xx(obj.uncertainIdx,:,:) = df_xx;
           dg_xx = squeeze(sum(dg_xx, 1));
           % Add the regularizer
           dg_xx = dg_xx + obj.Reg;
           %% Mixed partial derivative
           if nargout == 2
               % Deterministic terms
               nY = size(dP, 3);
               dz_y = sum(dP.*x',2);
               dz_y = dz_y + reshape(dw, [nX,1,nY]);
               dg_xy = 2 * (dz_y .* P + z .* dP);
               dg_xy(x < z,:,:) = 0;
               % Uncertain terms
               % Reshape to make the multiplications work
               dsigma_y = reshape(dsigma_y,size(dsigma_y,1),1,size(dsigma_y,2));
               dmu_y = reshape(dmu_y, size(dmu_y,1),1,size(dmu_y,2));
               dp_y = reshape(dp_y, size(dp_y, 1), 1, size(dp_y, 2));
               dc_y = reshape(dc_y, size(dc_y, 1), 1, size(dc_y, 2));
               % Calculate the terms that multiply each of the distribution
               % values for the mixed partials
               pdf_mult_y = 2.*(dsigma_y .* dsigma_x + m_sigma.*dsigma_xy).*(x_s + m_mu) + ...
                   2.*m_sigma.*dsigma_x.*dmu_y + 2.*m_sigma.*dsigma_y.*(del_nj + dmu_x) + m_sigma.^2 .* dmu_xy;
               dp_x_mult_y = 2.*m_sigma.*dsigma_y .*(x_s + m_mu) + m_sigma.^2 .* dmu_y;
               cdf_mult_y = 2.*(dsigma_y .* dsigma_x + m_sigma.*dsigma_xy + dmu_y .* dmu_x + m_mu.*dmu_xy);
               dc_x_mult_y = 2.*(m_sigma.*dsigma_y + m_mu .* dmu_y);
               % Calculate the mixed partial
               df_xy = -pdf_mult_y .* pdf - pdf_common .* dp_y - dp_x_mult_y .* dp_x - ...
                   dp_x_common .* dp_xy + cdf_mult_y .* cdf + cdf_common.* dc_y + dc_x_mult_y .* dc_x + ...
                   dc_x_common .* dc_xy;
               % Correct for degeneracy
               %df_xy(zeroIdx,:,:) = 2.*(dmu_y(zeroIdx,:,:) .* dmu_x(zeroIdx,:) + m_mu(zeroIdx) .* dmu_xy(zeroIdx,:,:));
               % Combine all the partial derivatives together to get the
               % derivatives of the cost
               dg_xy(obj.uncertainIdx,:,:) = df_xy;
               dg_xy = squeeze(sum(dg_xy, 1));
           end
       end
       function [pdf, cdf, dp_x, dc_x, dp_xx, dc_xx, dp_y, dc_y, dp_xy, dc_xy] = evalDistribution(obj, x, P, w, dP, dw)
           %% EVALDISTRIBUTION: Evaluate the Gaussian Distribution and its Derivatives
           %
           %    evalDistribution evaulates the Gaussian Distribution
           %    underlying the ERM problem and returns the value of the
           %    PDF, the value of the CDF, and various derivatives of those
           %    values. 
           %
           %    evalDistribution calculates the number of decision
           %    variables N, the number of stochastic decision variables
           %    Ns, and the number of additional parameters M internally
           %    from the arguments.
           %
           %   Arguments:
           %       OBJ:    an instance of the GaussianERM class
           %       x:      Nx1 double, decision variables of the ERM
           %               problem
           %       P:      NxN double, the LCP problem matrix
           %       w:      Nx1 double, the LCP problem vector
           %       dP:     NxNxM double, array of derivatives of P
           %       dw:     NxM double, array of derivatives of w
           %
           %   Return Values:
           %        pdf:    Ns x 1 double, the normal distribution evaluated
           %                at the stochastic decision variables
           %        cdf:    Ns x 1 double, the cdf of the normal
           %               distribution evaluated at the stochastic decision
           %               variables
           %        dp_x:   Ns x N, the derivative of pdf with respect to x
           %        dc_x:   Ns x N, the derivative of cdf with respect to x
           %        dp_y:   Ns x M, the derivative of pdf with respect to
           %                the addition parameters y
           %        dc_y:   Ns x M, the derivative of cdf with respect to
           %                additional parameters y
           %        dp_xx:  Ns x N x N, the second derivatives of pdf with
           %                respect to x
           %        dc_xx:  Ns x N x N, the second derivatives of cdf with
           %                respect to x
           %        dp_yx:  Ns x N x M, the mixed second derivatives of pdf
           %                with respect to x and y
           %        dc_yx:  Ns x N x M, the mixed second derivatives of cdf
           %                with respect to x and y
           
           % Get the mean and standard deviation for each problem
           if nargout > 6
               [m_mu, dmu_x, dmu_xx,dmu_y, dmu_xy] = obj.ermMean(x, P, w, dP, dw);
               [m_sigma, dsigma_x, dsigma_xx,dsigma_y, dsigma_xy] = obj.ermDeviation(x, P, w, dP, dw);
           else
               [m_mu, dmu_x, dmu_xx] = obj.ermMean(x, P, w);
               [m_sigma, dsigma_x, dsigma_xx] = obj.ermDeviation(x, P, w);
           end
           % Check for degeneracy
           zeroIdx = (m_sigma == 0);
           
           % Re-calculate the cost for the stochastic variables
           x_s = x(obj.uncertainIdx,:);

           % Normalize the decision variable
           tau = (x_s - m_mu)./m_sigma;
           % Calculate the PDF and CDF values
           pdf = normpdf(x_s, m_mu, m_sigma);
           cdf = normcdf(x_s, m_mu, m_sigma);
           % Degenerate case
           pdf(zeroIdx) = 0;
           cdf(zeroIdx) = 1;
           
           if nargout > 2
               % Calculate the first derivatives
               tau_const = tau.*dsigma_x + dmu_x;   %A common term in the derivatives of TAU
               tau_const(:,obj.uncertainIdx) = tau_const(:,obj.uncertainIdx) - eye(sum(obj.uncertainIdx));
               % The derivative of TAU wrt x
               dtau_x = - tau_const ./ m_sigma;
               
               % The derivatives of the pdf and cdf values wrt x
               dp_const = dsigma_x ./m_sigma + tau.*dtau_x;  % Term common in derivatives of pdf
               dp_x = - pdf .* dp_const;
               dc_x = - pdf .* tau_const;
               
               % Correct for degeneracy
               dp_x(zeroIdx,:) = 0;
               dc_x(zeroIdx,:) = 0;
               
               if nargout > 4
                   
                   % Calculate the second derivatives
                   % Broadcast some of the variables so the array sizes are
                   % consistent after multiplication
                   
                   dsigma_xk = reshape(dsigma_x, [size(dsigma_x, 1), 1, size(dsigma_x, 2)]);
                   dtau_xk = reshape(dtau_x, [size(dtau_x, 1), 1, size(dtau_x, 2)]);
                   % The second derivative of TAU
                   dtau_xx = tau_const .* dsigma_xk .* (m_sigma.^-1).^2 - (dtau_xk .* dsigma_x + tau.*dsigma_xx + dmu_xx)./m_sigma;
                   
                   dp_x2 = reshape(dp_x, [size(dp_x,1), 1, size(dp_x,2)]);

                   % Calculate the second derivatives
                   dp_xx = -dp_x2 .* dp_const - pdf .* (-dsigma_x .* dsigma_xk .* (m_sigma.^-1).^2 + dsigma_xx./m_sigma + dtau_x .* dtau_xk + tau.*dtau_xx);
                   dc_xx = dsigma_xk .* pdf .* dtau_x + m_sigma .* (dp_x2 .* dtau_x + pdf .* dtau_xx);
                   % Correct for degeneracy
                   dp_xx(zeroIdx,:,:) = 0;
                   dc_xx(zeroIdx,:,:) = 0;
                   if nargout > 6
                       % Broadcast some of the variables
                       dsigma_y = reshape(dsigma_y, [size(dsigma_y, 1), 1, size(dsigma_y, 2)]);
                       dmu_y = reshape(dmu_y, [size(dmu_y, 1), 1, size(dmu_y, 2)]);
                       % Calculate parameter derivatives
                       dtau_y = -(tau .* dsigma_y + dmu_y)./m_sigma;
                       dp_y = -pdf .* (dsigma_y ./ m_sigma + tau.*dtau_y);
                       dc_y = m_sigma .* pdf .* dtau_y;
                       % Correct for degeneracy
                       dp_y(zeroIdx,:,:) = 0;
                       dc_y(zeroIdx,:,:) = 0;
                       % the mixed derivative of TAU
                       dtau_xy = tau_const .* dsigma_y ./ (m_sigma.^2) - (dtau_y .* dsigma_x + tau.*dsigma_xy + dmu_xy)./m_sigma;
                       % The second derivatives of pdf and cdf
                       dp_xy = -dp_y .* dp_const - pdf .* (-dsigma_x .* dsigma_y ./ m_sigma.^2 + dsigma_xy ./m_sigma + dtau_x .* dtau_y + tau.*dtau_xy);
                       dc_xy = dsigma_x .* pdf .* dtau_y + m_sigma .* (dp_y .* dtau_x + pdf .* dtau_xy);
                       % Correct for degeneracy
                       dp_xy(zeroIdx,:,:) = 0;
                       dc_xy(zeroIdx,:,:) = 0;
                       % Squeeze down dp_y and dp_c
                       dp_y = reshape(dp_y, size(dp_y, 1), size(dp_y, 3));
                       dc_y = reshape(dc_y, size(dc_y, 1), size(dc_y, 3));
                   end
               end
           end
           
       end
       function [dg_x, dg_y] = costGradients(obj, x, P, w, dP, dw)
          %% COSTGRADIENTS: Returns the derivatives of the cost function 
          %
          % costGradients returns the values of the derivatives of the cost
          % function with respect to the decision variables x and
          % additional parameters y
          %
          %   Arguments:
          %       OBJ:    an instance of the GaussianERM class
          %       x:      Nx1 double, decision variables of the ERM
          %               problem
          %       P:      NxN double, the LCP problem matrix
          %       w:      Nx1 double, the LCP problem vector
          %       dP:     NxNxM double, array of derivatives of P
          %       dw:     NxM double, array of derivatives of w
          %
          %  Return values:
          %       dg_x:  1xN double, the derivative of the cost wrt the
          %                 decision variables
          %       dg_y:  1xM double, the derivative of the cost wrt the
          %                 additional parameters
          %
          % note: costGradients is normally used to verify the mixed
          % partial derivatives of the cost when changing the parameters is
          % not available (i.e. unit testing)
          
          % Calculate the NCP cost for the entire problem
          z = P * x + w;
          % Re-calculate the cost for the stochastic variables
          % Get the mean and standard deviation for each problem
          [m_mu, dmu_x, ~, dmu_y] = obj.ermMean(x, P, w, dP, dw);
          [m_sigma, dsigma_x, ~, dsigma_y] = obj.ermDeviation(x, P, w, dP, dw);
          % Check for degeneracy
          %zeroIdx = (m_sigma == 0);
          % Calculate the ERM residual
          [pdf, cdf, dp_x, dc_x,~, ~, dp_y, dc_y] = obj.evalDistribution(x, P, w, dP, dw);
          
          % Calculate the derivative wrt the decision variable x
          nX = numel(x);
          % First derivative for deterministic variables
          dg_x = 2 * z .* P;
          Dx = diag(x);
          dg_x(x < z,:) = 2.*Dx(x < z,:);
          % First derivative for stochastic variables
          % Get the stochastic variables
          % Get the uncertain variables
          x_s = x(obj.uncertainIdx,:);
          del_nj = zeros(length(x_s),nX);
          del_nj(:,obj.uncertainIdx) = eye(sum(obj.uncertainIdx));
          % Calculate the terms multiplying the distribution values
          pdf_common = 2.*m_sigma.*dsigma_x .* (x_s + m_mu) + m_sigma.^2 .*(del_nj + dmu_x);
          dp_x_common = m_sigma.^2 .*(x_s + m_mu);
          cdf_common = 2*(m_sigma.*dsigma_x + m_mu.*dmu_x - x_s.*del_nj);
          dc_x_common = m_sigma.^2 + m_mu.^2 - x_s.^2;
          % The derivative of the stochastic cost
          df_x = 2.*Dx(obj.uncertainIdx,:) - pdf_common .* pdf - dp_x_common .* dp_x + ...
              cdf_common .* cdf + dc_x_common .* dc_x;
          % Correct for degeneracy
          %df_x(zeroIdx,:) = 2.*m_mu(zeroIdx).*dmu_x(zeroIdx,:);
          % Combine the two costs
          dg_x(obj.uncertainIdx,:) = df_x;
          dg_x = sum(dg_x, 1);
          % Add the regularizer
          dg_x = dg_x + x'*obj.Reg;
          
          % Calculate the derivative wrt the paramters y
          % Derivative of deterministic variables
          dz_y = sum(dP .* x', 2);
          dz_y = reshape(dz_y, [size(dz_y, 1), size(dz_y, 3)]) + dw;
          dg_y = 2 .* z .* dz_y;
          dg_y(x < z, :) = 0;
          % Derivative of stochastic variables
          p_mult = 2.*m_sigma.*dsigma_y .* (x_s + m_mu) + m_sigma.^2 .* dmu_y;
          c_mult = 2.*(m_sigma.*dsigma_y + m_mu.*dmu_y);
          df_y = - p_mult .* pdf - dp_x_common .* dp_y + c_mult .* cdf + dc_x_common .* dc_y;
          % Correct for degeneracy
          %df_y(zeroIdx,:) = 2.*m_mu.*dmu_y;
          % Combine the derivatives
          dg_y(obj.uncertainIdx,:) = df_y;
          dg_y = sum(dg_y, 1);
       end
   end
end