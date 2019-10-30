classdef LagrangianModel < handle
    
    methods (Abstract)       
        M = inertiaMatrix(self,q)
        C = coriolisMatrix(self,q,qdot)
        N = gravityMatrix(self,q)
    end
    methods
        function xdot = dynamics(self, x)
            % Get the configuration and generalized velocity
            L = length(x);
            q = x(:,1:L/2);
            qdot = x(:,L/2+1:end);
            % Get the system properties
            M = self.inertiaMatrix(q);
            C = coriolisMatrix(self,q);
            N = gravityMatrix(self,q);
            % Calculate the acceleration
            qddot = -M\(C*qdot + N);
            xdot = [qdot;qddot];
        end
    end
end

