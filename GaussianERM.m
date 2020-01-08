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
    
    properties
        mu;              %Mean of the distribution over the uncertain parameter (NOT the mean used by the final ERM formula)
        sigma;           %Standard deviation of the distribution over the uncertain parameters
        options;         %Options structure for FMINUNC
    end
    properties (SetAccess = protected)
        uncertainIdx;    %Logical index indicating which variables are affected by the uncertainty
    end
   methods (Abstract)
       [m_mu, dmu_f, dmu_y, dmu_ff, dmu_yf] = ermMean(obj, f, P, w, dP, dw);
       [m_sigma, dsigma_f, dsigma_y, dsigma_ff, dsigma_fy] = ermDeviation(obj, P, w, dP, dw);
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
           obj.options = optimoptions('fmincon','display','none','SpecifyObjectiveGradient',true);
       end
       function [f, r] = solve(obj, P, w, dP, dw)
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
           
           % Cost function
           cost = @(x) ermCost(obj, x, P, w, dP, dw);
           % Initial guess
           f0 = zeros(size(w));
           % Solve the problem
           [f, r, exitflag] = fmincon(cost, f0, -P, w, [], [], zeros(size(f0)), [], [], obj.options);
           if exitflag<=0
              warning('FMINCON might have terminated before a solution was found'); 
           end
       end
       function df = gradient(obj, f, P, w, dP, dw)
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
           [g_ff, g_fy] = obj.ermSecondDerivatives(f, P, w, dP, dw);
           % Solve the system of equations to get the gradient
           df = - g_ff\g_fy;
           
           % NOTE: Check if we need to zero out the components of the
           % gradient corresponding to the zero solution.
       end
       function [r, dr] = ermCost(obj, x, P, w, dP, dw)
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
           %       dP:     NxNxM double, array of derivatives of P
           %       dw:     NxM double, array of derivatives of w
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
           % Re-calculate the cost for the stochastic variables
           x_s = x(obj.uncertainIdx,:);
           % Get the mean and standard deviation for each problem
           [m_mu, dmu_x] = obj.ermMean(x, P, w, dP, dw);
           [m_sigma, dsigma_x] = obj.ermDeviation(x, P, w, dP, dw);
           % Calculate the ERM residual
           [pdf, cdf, dp_x, dc_x] = obj.evalDistribution(x,P,w,dP,dw);
           % Calculate the residuals for only the stochastic part
           r(obj.uncertainIdx,:) = x_s.^2 - m_sigma.^2 .* (x_s + m_mu).*pdf + ...
               (m_sigma.^2 + m_mu.^2 - x_s^2).*cdf;
           % Sum all the residuals together
           r = sum(r);
           
           %% Calculate the gradient
           if nargout == 2
               % Gradient for the deterministic variables
               dr_1 = 2 * x';
               dr_1(x > z) = 0;
               dr_2 = 2 * (z.*P); 
               dr_2 = sum(dr_2,1);
               dr_2(x < z) = 0;
               % Gradient for the stochastic variables
               dx_x = zeros(sum(obj.uncertainIdx), length(x));
               dx_x(:,obj.uncertainIdx) = 1;
               % Gradient
               ds = 2*x_s - (2 * m_sigma.* ((x_s + m_mu) .* dsigma_x) + m_sigma.^2.*dx_x).*pdf -...
                   m_sigma.^2.*(x_s + m_mu).*dp_x + (2*m_sigma .* dsigma_x + 2*m_mu .* dmu_x - 2*x_s).*cdf +...
                   (m_sigma.^2 + m_mu.^2 - x_s.^2).*dc_x;
               % Sum the parts to get the overall gradient
               dr = dr_1 + dr_2 + sum(ds, 1);
           end
       end
   end
   methods (Sealed)
       function [dg_xx, dg_yx] = ermSecondDerivatives(obj, x, P, w, dP, dw)
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
           
           [nX,~,nY] = size(dP);
           
           %% Derivatives for determinsitic variables
           % Calculate the slack variables
           z = P*x + w;
           % Terms where z < x
           P2 = reshape(P,[nX,1,nX]);   % Broadcast P to make multiplication easier
           dg_xx = 2 * P .* P2;
           dz_y = sum(dP.*x',2);
           dz_y = dz_y + reshape(dw, [nX,1,nY]);
           dg_yx = 2 * (dz_y .* P + z .* dP);
           % Terms where x < z
           x_idx = x < z;
           I = eye(nX);
           delta_3 = I .* reshape(I,[nX,1,nX]);
           dg_xx(x_idx,:,:) = 2*delta_3(x_idx,:,:);
           dg_yx(x_idx,:,:) = 0;
           
           %% Derivatives for stochastic variables
           % Get the uncertain variables
           x_s = x(obj.uncertainIdx,:);
           del_nj = zeros(length(x_s),nX);
           del_nj(:,obj.uncertainIdx) = eye(sum(obj.uncertainIdx));
           % Derivatives of the distribution values
           [pdf, cdf, dp_x, dc_x, dp_y, dc_y, dp_xx, dc_xx, dp_yx, dc_yx] = obj.evalDistribution(x, P, w, dP, dw);
           [m_sigma, dsigma_x, dsigma_y, dsigma_xx, dsigma_yx]  = obj.ermDeviation(x, P, w, dP, dw);
           [m_mu, dmu_x, dmu_y, dmu_xx, dmu_yx] = obj.ermMean(x, P, w, dP, dw);
           % Reshape to make the multiplications work
           dsigma_y = reshape(dsigma_y,size(dsigma_y,1),1,size(dsigma_y,2));
           dsigma_xk = reshape(dsigma_x,size(dsigma_x,1),1,size(dsigma_x,2));
           del_nk = reshape(del_nj,size(del_nj,1),1,size(del_nj,2));
           dmu_y = reshape(dmu_y, size(dmu_y,1),1,size(dmu_y,2));
           dmu_xk = reshape(dmu_x,size(dmu_x,1),1,size(dmu_x,2));
           dp_y = reshape(dp_y,size(dp_y,1),1,size(dp_y,2));
           dc_y = reshape(dc_y,size(dc_y,1),1,size(dc_y,2));
           dp_xk = reshape(dp_x,size(dp_x,1),1,size(dp_x,2));
           dc_xk = reshape(dc_x,size(dc_x,1),1,size(dc_x,2));
           % Calculate terms common to both derivatives
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
           % Calculate the terms that multiply each of the distribution
           % values for the mixed partials
           pdf_mult_y = 2.*(dsigma_y .* dsigma_x + m_sigma.*dsigma_yx).*(x_s + m_mu) + ...
               2.*m_sigma.*dsigma_x.*dmu_y + 2.*m_sigma.*dsigma_y.*(del_nj + dmu_x) + m_sigma.^2 .* dmu_yx;
           dp_x_mult_y = 2.*m_sigma.*dsigma_y .*(x_s + m_mu) + m_sigma.^2 .* dmu_y;
           cdf_mult_y = 2.*(dsigma_y .* dsigma_x + m_sigma.*dsigma_yx + dmu_y .* dmu_x + m_mu.*dmu_yx);
           dc_x_mult_y = 2.*(m_sigma.*dsigma_y + m_mu .* dmu_y);
           % Calculate the mixed partial
           df_yx = -pdf_mult_y .* pdf - pdf_common .* dp_y - dp_x_mult_y .* dp_x - ...
               dp_x_common .* dp_yx + cdf_mult_y .* cdf + dc_x_common .* dc_y + dc_x_mult_y .* dc_x + ...
               dc_x_common .* dc_yx;
           %% Combine the derivatives of the deterministic and stochastic variables
           % Combine all the partial derivatives together to get the
           % derivatives of the cost
           dg_xx(obj.uncertainIdx,:,:) = df_xx;
           dg_yx(obj.uncertainIdx,:,:) = df_yx;
           dg_xx = squeeze(sum(dg_xx, 1));
           dg_yx = squeeze(sum(dg_yx, 1));
       end
       function [pdf, cdf, dp_x, dc_x, dp_y, dc_y, dp_xx, dc_xx, dp_yx, dc_yx] = evalDistribution(obj, x, P, w, dP, dw)
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
           
           % Re-calculate the cost for the stochastic variables
           x_s = x(obj.uncertainIdx,:);
           % Get the mean and standard deviation for each problem
           [m_mu, dmu_x, dmu_y, dmu_xx, dmu_xy] = obj.ermMean(x, P, w, dP, dw);
           [m_sigma, dsigma_x, dsigma_y, dsigma_xx, dsigma_xy] = obj.ermDeviation(x, P, w, dP, dw);
           % Normalize the decision variable
           tau = (x_s - m_mu)./m_sigma;
           % Calculate the PDF and CDF values
           pdf = normpdf(x_s, m_mu, m_sigma);
           cdf = normcdf(x_s, m_mu, m_sigma);
           
           if nargout > 2
               % Calculate the first derivatives
               tau_const = tau.*dsigma_x + dmu_x;   %A common term in the derivatives of TAU
               tau_const(:,obj.uncertainIdx) = tau_const(:,obj.uncertainIdx) - eye(sum(obj.uncertainIdx));
                % The derivative of TAU wrt x
               dtau_x = - tau_const ./ m_sigma;

               % The derivatives of the pdf and cdf values wrt x
               dp_const = dsigma_x ./m_sigma + dtau_x;  % Term common in derivatives of pdf
               dp_x = - pdf .* dp_const;
               dc_x = m_sigma .* pdf .* dtau_x;
               
               if nargout > 4
                   % Calculate the derivatives wrt the parameters y
                   dtau_y = (tau .* dsigma_y + dmu_y)./m_sigma;
                   dp_y = -pdf .* (dsigma_y ./ m_sigma + dtau_y);
                   dc_y = m_sigma .* pdf .* dtau_y;
                   
                   % Calculate the second derivatives
                   % Broadcast some of the variables so the array sizes are
                   % consistent after multiplication
                   dsigma_y = reshape(dsigma_y, [size(dsigma_y, 1), 1, size(dsigma_y, 2)]);
                   dsigma_x2 = reshape(dsigma_x, [size(dsigma_x, 1), 1, size(dsigma_x, 2)]);
                   dtau_y = reshape(dtau_y, [size(dtau_y, 1), 1, size(dtau_y, 2)]);
                   dp_x2 = reshape(dp_x, [size(dp_x,1), 1, size(dp_x,2)]);
                   dp_y2 = reshape(dp_y, [size(dp_y,1), 1, size(dp_y,2)]);
                   
                   % the second derivatives of TAU 
                   dtau_xx = tau_const .* dsigma_x2 ./ (m_sigma.^2) - (dtau_x .* dsigma_x2 + tau.*dsigma_xx + dmu_xx)./m_sigma;
                   dtau_yx = tau_const .* dsigma_y ./ (m_sigma.^2) - (dtau_y .* dsigma_x + tau.*dsigma_xy + dmu_xy)./m_sigma;
                   
                   % The second derivatives of pdf and cdf
                   dp_xx = -dp_x2 .* dp_const - pdf .* (dsigma_x .* dsigma_x2 ./ m_sigma.^2 + dsigma_xx./m_sigma + dtau_xx);
                   dp_yx = -dp_y2 .* dp_const - pdf .* (dsigma_x .* dsigma_y ./ m_sigma.^2 + dsigma_xy + dtau_yx);
                   
                   dc_xx = dsigma_x2 .* pdf .* dtau_x + m_sigma .* (dp_x2 .* dtau_x + pdf .* dtau_xx);
                   dc_yx = dsigma_x .* pdf .* dtau_y + m_sigma .* (dp_y2 .* dtau_x + pdf .* dtau_yx);
               end
           end
           
       end
   end
end