classdef logisticREM < ContactSolver
    
   properties
       opts;
       sigma;
   end
   methods
       function obj = logisticREM(sig)
          obj.sigma = sig;
          obj.opts = optimoptions('lsqnonlin','MaxFunctionEvaluations',2000,'StepTolerance',1e-8);
       end
       function [x, r] = solve(obj, P, w)
           residuals = @(x) obj.remResiduals(x, P, w);
           x0 = zeros(size(w));
           x = lsqnonlin(residuals,x0,[],[],obj.opts);
           r = obj.remCost(x, P, w);
       end
       function g = remCost(obj, x, P, w)
           r = obj.remResiduals(x, P, w);
           g = sum(r.^2);
       end
       function r = remResiduals(obj, x, P, w)
           % SLCP Residuals with Logistic function
           z = P * x + w;
           p = exp((x - z)./obj.sigma);
           r = a - obj.sigma .* log(1+p);
       end
       function df = gradient(obj, f, P, z, dP, dz)
           % Gradient: Gradient of the solution wrt other problem
           % variables
           %   TODO: Implement gradient for NCP solver
           warning("Gradient for logisticREM hasn't been implemented yet");
           df = zeros(length(f), size(dP,3));
       end
   end
   
end