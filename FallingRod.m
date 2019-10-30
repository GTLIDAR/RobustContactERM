classdef FallingRod < ContactDynamics
   
    
    properties
        mass = 1;
        length = 1;
        inertia = 1;
        radius = 1;
        friction = 1;
        timestep = 1;
    end
    methods
      
        function self = FallingRod(m, l, r, J, mu, h)
           self.mass = m;
           self.length = l;
           self.radius = r;
           self.inertia = J;
           self.friction = mu;
           self.timestep = h;
        end
        function M = inertiaMatrix(self, ~)
            M = zeros(3,3);
            M(1,1) = self.mass;
            M(2,2) = self.mass;
            M(3,3) = self.inertia;
        end
        function C = coriolisMatrix(~,~,~)
            C = zeros(3,3);
        end
        function N = gravityMatrix(self,~)
            N = zeros(3,1);
            N(2) = self.mass * 9.81;
        end
        function [normals, alphas] = contactNormal(self, q)
            normals = zeros(2,3);
            normals(1,:) = [0, 1, -self.length/2 * cos(q(3))];
            normals(2,:) = [0, 1,  self.length/2 * cos(q(3))];
            
            alphas = zeros(2,1);
            alphas(1,:) = self.radius + self.length/2 * (sin(q(3)) - q(3)*cos(q(3)));
            alphas(2,:) = self.radius + self.length/2 * (q(3)*cos(q(3)) - sin(q(3)));        
        end
        function [Jn, Jt] = contactJacobian(self, q)
            
            Jt = zeros(3, 4);
            Jn = zeros(3, 2);
            
            Jt(1,:) = [1,-1,1,-1];
            Jt(3,:) = self.length/2 * [sin(q(3)), -sin(q(3)), -sin(q(3)), sin(q(3))];
            Jn(2,:) = [1,1];
            Jn(3,:) = self.length/2 * [-cos(q(3)), cos(q(3))];
        end
        function [x,y] = cartesian(self, q)
            x = q(1) + (self.length/2  + self.radius) * cos(q(3)) * [1,-1];
            y = q(2) + (self.length/2  + self.radius) * sin(q(3)) * [1,-1];
        end
        function plotTrajectory(self,t, x, f)
            % Velocity figure
            figure();
            labels = {'Horizontal','Vertical','Angular'};
            for n = 1:3
                subplot(3,1,n);
                plot(t, x(n+3, :));
                ylabel(labels{n});
            end
            xlabel('Time (s)');
            subplot(3,1,1);
            title('Velocities');
            % Force figure
            % Re-map the forces
            A = [1, -1, 0, 0, 0, 0;
                0, 0, 1, -1, 0, 0;
                0, 0, 0, 0,  1, 0;
                0, 0, 0, 0,  0, 1];
            f = A*f;
            % Plot the forces
            labels = {'Friction 1','Friction 2','Normal 1','Normal 2'};
            figure();
            for n = 1:4
                subplot(4,1,n);
                plot(t,f(n,:));
                ylabel(labels{n});
            end
            xlabel('Time (s)');
            subplot(4,1,1)
            title('Contact Forces');
            % Plot the trajectory of the rod as it falls
            figure();
            for n = 1:numel(t)
                [xp,yp] = self.cartesian(x(:,n));
                plot(xp,yp);
                hold on;
            end
        end
        function frames = animate(self,x)
           figure();
           q = x(1:3,:);
           ax = gca;
           [x0,y0] = self.getPoints([0,0,0]');
           drawing = plot(ax,x0,y0,'b-','LineWidth',1.5);
           hold on;
           plot([-1,1],[0,0],'k-','LineWidth',2);
           xlim([-1,1]);
           ylim([-0.5,1.5])
           frames(1:size(x,2)) = getframe();
           for n = 1:size(q,2)
               [x,y] = self.getPoints(q(:,n));
               set(drawing,'XData',x,'YData',y);
               drawnow;
               frames(n) = getframe();
           end
        end
    end
    methods (Access = private)
        function X = drawHome(self)
            % Return homogeneous coordinates of the points on the body (for
            % drawing the rigid body)
            
            % Rod starts with center of mass at (x,y) = (0,0), and with
            % Orientation theta = 0
            
            % Centers of the semicircle ends
            LCenter = [-self.length/2, 0];
            RCenter = [self.length/2, 0];
            % Create the circular ends
            L_Th = linspace(pi/2, 3*pi/2, 20);
            R_Th = linspace(-pi/2, pi/2, 20);
            
            LCircle = zeros(2,20);
            LCircle(1,:) = LCenter(1) + self.radius * cos(L_Th);
            LCircle(2,:) = LCenter(2) + self.radius * sin(L_Th);
            
            RCircle = zeros(2,20);
            RCircle(1,:) = RCenter(1) + self.radius * cos(R_Th);
            RCircle(2,:) = RCenter(2) + self.radius * sin(R_Th);
            
            X = [LCircle, RCircle, LCircle(:,1)];
            % Append X with 1s for homogeneous coordinates
            X = [X; ones(1,size(X,2))];
            
        end
        function g = transform(~,q)
            % Create the rigid body transform from the home configuration
            % to the current configuration
            g = eye(3);
            g(1:2,1:2) = [cos(q(3)), -sin(q(3)); 
                          sin(q(3)), cos(q(3))];
            g(1:2,3) = q(1:2);
        end
        function [x,y] = getPoints(self,q)
            % Get the points in a standard configuration, in homogeneous
            % coordinates
            X = self.drawHome();
            % Apply the transformation specified by the current
            % configuration
            g = self.transform(q);
            X = g*X;
            % Extract the cartesian points
            x = X(1,:);
            y = X(2,:);
        end
    end
end