classdef PathLCPSolver < ContactSolver
    
    methods (Static)
        function [f,r] = solve(P,z)
            try
                f = pathlcp(P,z);
                r = f'*(P*f + z);
            catch
                f = zeros(size(z));
                r = f'*(P*f + z);
                warning('PATH Solver failed');
            end
        end
        function  df = gradient(f, P, ~, dP, dz)
            % Get only the parts of the problem corresponding to nonzero
            % solutions
            df = zeros(size(dz));
            id = ~(f==0);
            if ~isempty(id) 
                dz = dz(id,:);
                P = P(id, id);
                dP = dP(id,id,:);
                f = f(id);
                % Tensor multiplication
                dP_f = times(dP, reshape(f, [1, numel(f),1]));
                dP_f = sum(dP_f,2);
                dP_f = reshape(dP_f, size(dP_f, 1), size(dP_f, 3));
                % Nonzero gradients
                df(id,:) = - pinv(P) * (dP_f + dz);
            end
        end
    end
end