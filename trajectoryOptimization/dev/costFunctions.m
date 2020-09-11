classdef costFunctions

    %% COSTFUNCTIONS: a suite of parameterized cost functions for use in optimal control problems
    %   costFunctions is a static class containing several common, mostly
    %   quadratic cost functions that are compatible with Drake's
    %   trajectory optimization tools. 
    %
    %   Many of the cost functions presented are parameterized by a
    %   reference value and a weighting matrix
    
    
    %   Luke Drnach
    %   July 31, 2020
    
    methods (Static)
        function [g,dg] = quadraticControl(h,x,u,ref,W)
            %% Quadratic cost penalty on control effort
            % Deviation from reference
            du = u -ref;
            % Weighted Cost
            g = 0.5 * du'*W*du;
            % Differential cost
            dg = [zeros(1,numel(h) + numel(x)), du'*W];
        end
        
        function [g,dg] = quadraticState(h,x,u,ref,W)
            %% QuadraticState: Quadratic cost on state deviation
            %
            %   Penalizes deviation from reference state
            
            % Deviation from reference
            dx = x - ref;
            % Weighted cost
            g = 0.5 * dx'*W*dx;
            % Differential cost
            dg = [zeros(1,numel(h)), dx'*W, zeros(1,numel(u))];
        end
  
        function [g,dg] = quadraticForce(h,l,ref,W)
           %% 
           dl = l - ref;
           g = 0.5*dl'*W*dl;
           dg = [zeros(1,numel(h)), dl'*W];
        end
        function [g,dg] = minimumTime(h,x)
            
            g = sum(h);
            dg = [ones(1,numel(h)), zeros(1,numel(x))];
            
        end
        function [g,dg] = terminalState(h,x,ref,W)
            
            dx = x - ref;
            g = 0.5 * dx' * W * dx;
            dg = [zeros(1,numel(h)), dx'*W];
            
        end
        function [g,dg] = smoothState(dt,dx,du,dl,W)
            
            g = 0.5 * dx'*W*dx;
            dg =  [zeros(1,numel(dt)), dx'*W, zeros(1,numel(du) + numel(dl))];
            
        end
        function [g,dg] = smoothControl(dt, dx, du, dl, W)
            
            g = 0.5 * du'*W*du;
            dg = [zeros(1,numel(dt)+numel(dx)), du'*W, zeros(1,numel(dl))];
            
        end
        function [g,dg] = smoothForce(dt,dx,du,dl,W)
            
            g = 0.5 * dl'*W*dl;
            dg = [zeros(1,numel(dt)+numel(du) + numel(dx)), dl'*W];
            
        end
    end
end