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
            
            % Tensor multiplication
%             dP_f = times(dP, reshape(f, [1, numel(f), 1]));
%             dP_f = squeeze(sum(dP_f, 2));
%             % Calculate gradient
%             df = - P \ (dP_f + dz);
%             % When f = 0, df must also be 0
%             df(f==0, :) = 0;
              % Get only the nonzero indices
              df = zeros(size(dz));
              id = ~(f==0);
              if ~isempty(id)
                  dz = dz(id,:);
                  P = P(id, id);
                  dP = dP(id,id,:);
                  % Tensor multiplication
                  f = f(id);
                  dP_f = times(dP, reshape(f, [1, numel(f),1]));
                  dP_f = squeeze(sum(dP_f,2));
                  % Calculate gradient
                  df(id,:) = - pinv(P) * (dP_f + dz);
              end
              
        end
    end
end