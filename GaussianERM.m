classdef GaussianERM < ContactSolver
    
   properties
       mu;
       sigma;
       uncertainIdx;
       options;
   end
   methods (Abstract)
       [m_mu, dmu] = ermMean(obj, P, w, dP, dw);
       [m_sigma, dsigma] = ermDeviation(obj, P, w, dP, dw);
   end
   
   methods
       function obj = GaussianERM(m, s, idx)
           obj.mu = m;
           obj.sigma = s;
           obj.uncertainIdx = idx;
           obj.options = optimoptions('fminunc','Algorithm','trust-region',...
               'display','none','SpecifyObjectiveGradient',true);
       end
       function [f, r] = solve(obj, P, w)
           
           
           % Cost function
           cost = @(x) ermCost(obj, x, P, w, dP, dw);
           % Initial guess
           f0 = zeros(size(w));
           % Solve the problem
           [f, r, exitflag] = fminunc(cost, f0, obj.options);
           if exitflag<=0
              warning('FMINUNC might have terminated before a solution was found'); 
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
           
           % Form the linear system describing the derivatives assuming
           % all variables are deterministic
           A = P;
           dP_f = times(dP, reshape(f,[1, numel(f), 1]));
           dP_f = sum(dP_f, 2);
           dP_f = reshape(dP_f, [size(dP_f,1),size(dP_f,3)]);
           b = -(dP_f + dw);
           % Recalculate the components relating to the stochastic
           % variables
           [As,bs] = obj.ermGradientCoefficients(f, P, w, dP, dw);
           A(obj.uncertainIdx,:) = As;
           b(obj.uncertainIdx,:) = bs;
           % Remove the components for which is the solution is zero
           nonzeroIdx = (f ~= 0);
           A = A(nonzeroIdx,nonzeroIdx);
           b = b(nonzeroIdx,:);
           % Solve the linear system for the derivatives
           df_nz = A\b;
           df = zeros(length(f), size(df_nz, 2));
           df(nonzeroIdx,:) = df_nz;
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
           
           pdf = normpdf(x_s, m_mu, m_sigma);
           cdf = normcdf(x_s, m_mu, m_sigma);
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
               % Derivatives of the pdf and cdf values
               tau = (x_k - m_mu)./m_sigma;
               dtau_x = - dsigma_x .* (tau./m_sigma) + (dx_x - dmu_x)./m_sigma;
               dp_x = -(pdf./m_sigma) .* (dsigma_x ./ m_sigma + dtau_x);
               dc_x = dtau_x .* (m_sigma .* pdf);
               % Gradient
               ds = 2*x_k - (2 * m_sigma.* ((x_k + m_mu) .* dsigma_x) + m_sigma.^2.*dx_x).*pdf -...
                   m_sigma.^2.*(x_k + m_mu).*dp_x + (2*m_sigma .* dsigma_x + 2*m_mu .* dmu_x - 2*x_k).*cdf +...
                   (m_sigma.^2 + m_mu.^2 - x_k.^2).*dc_x;
               % Sum the parts to get the overall gradient
               dr = dr_1 + dr_2 + sum(ds, 1);
           end
       end
   end
   methods (Sealed)
       function [A, b] = ermGradientCoefficients(obj, f, P, w, dP, dw)
           
           
           %% Some set-up
           % Get the solution values that correspond to the uncertainty.
           f_k = f(obj.uncertainIdx,:);
           % Get the parameters of the ERM problem
           [m_mu, dmu_f, dmu_y] = obj.ermMean(f, P, w, dP, dw);
           [m_sigma, dsigma_f, dsigma_y] = obj.ermDeviation(f,P,w,dP,dw);
           %% Calculate the probability values and their gradients
           tau = (f_k - m_mu)./m_sigma;
           df_f = zeros(sum(obj.uncertainIdx), length(f));
           df_f(:,obj.uncertainIdx) = 1;
           pdf = normpdf(f(obj.uncertainIdx), m_mu, m_sigma);
           cdf = normcdf(f(obj.uncertainIdx), m_mu, m_sigma);
           % Derivatives of the normalized values, tau
           dtau_f = - dsigma_f .* (tau./m_sigma) + (df_f - dmu_f)./m_sigma;
           dtau_y = - dsigma_y .*(tau./m_sigma) - dmu_y ./ m_sigma;
           % Derivatives of the probability density function, pdf
           dp_f = -(pdf./m_sigma) .* (dsigma_f ./ m_sigma + dtau_f);
           dp_y = -(pdf./m_sigma) .* (dsigma_y ./ m_sigma + dtau_y);
           % Derivatives of the cumulative density function, cdf
           dc_f = dtau_f .* (m_sigma .* pdf);
           dc_y = dtau_y .* (m_sigma .* pdf);
           %% Calculate the coefficients of the gradients
           % The coefficients of the gradients
           A = (dsigma_f .* (2 .*m_sigma.*(f_k - m_mu)) + m_sigma.^2.*dmu_f).*pdf -...
               m_sigma.^2 .* (f_k + m_mu) .* dp_f + (dsigma_f .* (2*m_sigma) + dmu_f .*(2.*m_mu)).*cdf + ...
               (m_sigma.^2 + m_mu.^2 - f_k.^2).*dc_f;
           A(:,obj.uncertainIdx) = A(:, obj.uncertainIdx) - (m_sigma.^2 .* pdf + ...
               2*f(obj.uncertainIdx,:) .* (1 - cdf));
           % The constant coefficients
           b = - (dsigma_y .*(2.*m_sigma).*(f_k + m_mu) + dmu_y .*(m_sigma.^2)).*pdf - ...
               dp_y .* m_sigma.^2 .* (f_k + m_mu) + (dsigma_y .* (2.*m_sigma) + dmu_y .* (2.*m_mu)).*cdf + ...
               dc_y .* (m_sigma.^2 + m_mu.^2 - f_k.^2);
       end
   end
end