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
            q = x(1:L/2);
            qdot = x(L/2+1:end);
            % Get the system properties
            M = self.inertiaMatrix(q);
            C = coriolisMatrix(self,q, qdot);
            N = gravityMatrix(self,q);
            % Calculate the acceleration
            qddot = -M\(C*qdot + N);
            xdot = [qdot;qddot];
        end
        function [x,t] = simulate(self,x0,T, dt)
            % Pre-initialize arrays
            nX = numel(x0);
            time = 0:dt:T;
            nT = numel(time);
            x= zeros(nX,nT);
            x(:,1) = x0;
            % Loop over the dynamics, using Euler integration to calculate
            % the next state
            for t = 1:nT-1
               xdot = self.dynamics(x(:,t));
               x(:,t+1) = x(:,t) + dt*xdot;
            end
        end
        function animate(self,x,dt)
            %
            %
            %   NOTE: Here, I use explicit Euler integration to simulate
            %   the model forward. This scheme is both numerically unstable
            %   and does not preserve the energy/momenta of the system.
            %   Simulations are likely to be very inaccurate, especially
            %   over longer time horizons.
            f = figure();
            ax = gca;
            nQ = numel(x)/2;
            % Draw the ininitial figure
            self.draw(x(1:nQ), ax);
            % Run while the figure is open
            while(ishandle(f))
                xdot = self.dynamics(x);
                % Use a sympletic euler integration scheme
                %x(nQ+1:end) = xdot(nQ+1:end)*dt + x(nQ+1:end);
                %x(1:nQ) = x(1:nQ) + x(nQ+1:end)*dt;
                x = x + dt*xdot;   % Explicit Euler - Do NOT USE
                % Update the figure
                [xp,yp] = self.draw(x(1:nQ));
                set(get(ax,'Children'),'XData',xp,'YData',yp);
                drawnow;
            end
        end
    end
end

