classdef MINSolver < ContactSolver
    
    properties
       opts; 
    end
    
    methods
        function obj = MINSolver()
            obj.opts = optimoptions('lsqnonlin','MaxFunctionEvaluations',2000,'StepTolerance',1e-8);
        end
        function [x, r] = solve(obj, P, w)
            % Use the MIN function to solve the LCP (via lsqnonlin)
            residuals = @(x) obj.ncpResiduals(x, P, w);
            x0 = zeros(size(w));
            x = lsqnonlin(residuals,x0,[],[],obj.opts);
            r = obj.ncpCost(x, P, w);
        end
        function g = ncpCost(obj, x, P, w)
            % COST: Sum-of-squared residuals
            r = obj.ncpResiduals(x, P, w);
            g = sum(r.^2);
        end
        function r = ncpResiduals(obj, x, P, w)
            % NCP Residual function using MIN
            z = P*x + w;
            r = min(x, z);
        end
        function df = gradient(obj, f, P, z, dP, dz)
            % Gradient: Gradient of the solution wrt other problem
            % variables
            %   TODO: Implement gradient for NCP solver
            warning("Gradient for MINSolver hasn't been implemented yet");
            df = zeros(length(f), size(dP,3));
        end
    end
end
