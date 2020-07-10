classdef FootedHopper <  Manipulator & DifferentiableContactDynamics 
    
    properties
        baseMass(1,1) double {mustBePositive, mustBeFinite} = 10;
        masses(3,1) double {mustBePositive, mustBeFinite} = [1,1,1/10];
        lengths(3,1) double {mustBeNonnegative, mustBeFinite} = [1, 1, 1/10];
        com(3,1) double = [1/2, 1/2, 1/2];
        inertias(3,1) double = [1/12, 1/12, 1/(12 * 10^3)];
        heelFraction(1,1) double = 1/2;
    end
    properties (Hidden)
        base_radius(1,1) double = 0.25;
        link_radius(1,1) double = 0.02;
        numQ = 5;
        numU = 3;
    end
    
methods
    function self = FootedHopper(varargin)
        %% SingleLegHopper: Class Constructor
        
        nQ = 5;
        nU = 3;
        self = self@Manipulator(nQ, nU);
                
        if nargin > 0
            if ~isempty(varargin{1})
            	self.masses = varargin{1};
            end
            if nargin > 1
                if ~isempty(varargin{2})
                    self.lengths = varargin{2};
                end
                if nargin > 2
                    if ~isempty(varargin{3})
                        self.com = varargin{3};
                    end
                    if nargin > 3
                        if ~isempty(varargin{4})
                            self.inertias = varargin{4};
                        end
                    end
                end
            end
        end
    end
    %% --------------- KINEMATIC PROPERTIES -------------------  %%
    function [x,y] = positions(self,q)
       % Returns the positions of the links in the system

              
       x = [q(1), self.lengths(1) * sin(q(3)), self.lengths(2)*sin(q(3)+q(4)),0,0];
       x = cumsum(x);
       x(4) = x(4) +  (1 - self.heelFraction)*self.lengths(3)*sin(q(3)+q(4)+q(5));
       x(5) = x(5) - self.heelFraction*self.lengths(3)*sin(q(3)+q(4)+q(5));
       
       y = [q(2), -self.lengths(1) * cos(q(3)), -self.lengths(2)*cos(q(3) + q(4)), 0, 0];
       y = cumsum(y);
       y(4) = y(4) - (1 - self.heelFraction)*self.lengths(3)*cos(q(3)+q(4)+q(5));
       y(5) = y(5) + self.heelFraction*self.lengths(3)*cos(q(3)+q(4)+q(5));

       
       
    end
    function x = kinematics(self, q)
       % Returns the position of the endpoint

       [x,y] = self.positions(q);
       x = [x(4:5); y(4:5)];
    end
    function [x, dx] = footPosition(self,q)
        % Returns the position of the joint connecting the foot to the leg
        
        % Initialize
        x = zeros(2,1);
        dx = zeros(2,numel(q));
        % Positions
        x(1) = q(1) + self.lengths(1)*sin(q(3)) + self.lengths(2)*sin(q(3)+q(4));
        x(2) = q(2) - self.lengths(1)*cos(q(3)) - self.lengths(2)*cos(q(3)+q(4));
        
        % Gradients
        dx(1,:) = [1, 0, self.lengths(1)*cos(q(3)) + self.lengths(2)*cos(q(3)+q(4)), self.lengths(2)*cos(q(3)+q(4)), 0];
        dx(2,:) = [0, 1, self.lengths(1)*sin(q(3)) + self.lengths(2)*sin(q(3)+q(4)), self.lengths(2)*sin(q(3)+q(4)), 0];
        
    end
    function [J, dJ] = jacobian(self, q)
       % Returns the endpoint Jacobian for the system
       % The endpoint Jacobian relates the generalized velocities to the
       % Cartesian coordinate velocities of the end of the second pendulum
       
       
       %% Relative Positions of all revolute joints
       x = [self.lengths(1)*sin(q(3)), self.lengths(2)*sin(q(3)+q(4)), (1-self.heelFraction)*self.lengths(3)*sin(q(3)+q(4)+q(5)),...
           -self.heelFraction*self.lengths(3)*sin(q(3)+q(4)+q(5))];
       y = [self.lengths(1)*cos(q(3)), self.lengths(2)*cos(q(3)+q(4)), (1-self.heelFraction)*self.lengths(3)*cos(q(3)+q(4)+q(5)),...
           -self.heelFraction*self.lengths(3)*cos(q(3)+q(4)+q(5))];
       
       %% Toe Jacobian
       J1 = zeros(2,5);
       % Base 
       J1(1:2,1:2) = eye(2); 
       % Joints
       J1(1,3:5) = y(1:3);
       J1(2,3:5) = x(1:3);
       J1(:,5:-1:3) = cumsum(J1(:,5:-1:3), 2);
       %% Heel Jacobian
       J2 = zeros(2,5);
       J2(1:2,1:2) = eye(2);
       % Base
       % Joints
       J2(1,3:5) = y([1,2,4]);
       J2(2,3:5) = x([1,2,4]);
       J2(:,5:-1:3) = cumsum(J2(:,5:-1:3),2);
       %% Concatenate the Jacobians together 
       J = [J1;J2];
       %% Jacobian Derivatives
       if nargout == 2 
          dJ1 = zeros(2,5,5);
          dJ2 = zeros(2,5,5);
          
          % Toe Jacobian Derivatives
          dJ1(1,3:5,5) = -x(3);
          dJ1(2,3:5,5) = y(3);
          dJ1(:,:,4) = dJ1(:,:,5);
          dJ1(1,3:4,4) = dJ1(1,3:4,4) - x(2);
          dJ1(2,3:4,4) = dJ1(2,3:4,4) + y(2);
          dJ1(:,:,3) = dJ1(:,:,4);
          dJ1(1,3,3) = dJ1(1,3,3)-x(1);
          dJ1(2,3,3) = dJ1(2,3,3)+ y(1);
          
          % Heel Jacobian Derivative
          dJ2(1,3:5,5) = -x(4);
          dJ2(2,3:5,5) = y(4);
          dJ2(:,:,4) = dJ2(:,:,5);
          dJ2(1,3:4,4) = dJ2(1,3:4,4) - x(2);
          dJ2(2,3:4,4) = dJ2(2,3:4,4) + y(2);
          dJ2(:,:,3) = dJ2(:,:,4);
          dJ2(1,3,3) = dJ2(1,3,3)-x(1);
          dJ2(2,3,3) = dJ2(2,3,3)+ y(1);
          
          % Concatenate
          dJ = [dJ1; dJ2];
       end
       
    end
    %% ---------------- DYNAMIC PROPERTIES --------------------- %%
    function [M, dM, HM] = massMatrix(self, q)
       %% Returns the Generalized Mass Matrix for the System
       %
       %
       
       %% Calculate Inertial constants
       % Constants
       k1 = self.baseMass + sum(self.masses);
       k2 = self.inertias(3) + self.masses(3) * (self.com(3) - self.heelFraction)^2 * self.lengths(3)^2;
       k3 = self.inertias(2) + (self.masses(2) * self.com(2)^2 + self.masses(3)) * self.lengths(2)^2;
       k4 = self.inertias(1) + (self.masses(1)  * self.com(1)^2 + self.masses(2) + self.masses(3)) * self.lengths(1)^2;
       % Sine / Cosine Multipliers
       g1 = (self.masses(1)*self.com(1) + self.masses(2) + self.masses(3))*self.lengths(1);
       g2 = (self.masses(2)*self.com(2) + self.masses(3))*self.lengths(2);
       g3 = self.masses(3)*(self.com(3) - self.heelFraction)*self.lengths(3);
       m1 = g3 * self.lengths(1);
       m2 = g3 * self.lengths(2);
       m3 = g2 * self.lengths(1);
       
       %% Sine and Cosine
       s1 = sin(q(3));
       s2 = sin(q(4));
       s3 = sin(q(5));
       c1 = cos(q(3));
       c2 = cos(q(4));
       c3 = cos(q(5));
       
       s12 = sin(q(3) + q(4));
       s23 = sin(q(4) + q(5));
       s123 = sin(q(3) + q(4) + q(5));
       
       c12 = cos(q(3) + q(4));
       c23 = cos(q(4) + q(5));
       c123 = cos(q(3) + q(4) + q(5));
       
       %% Construct the Mass Matrix
       % Diagonals
       dM = [k1, k1, k2 + k3 + k4 + 2*(m1 * c23 + m2 * c3 + m3 * c2), k2 + k3 + 2 * m2 * c3, k2];
       % Upper Triangle
       M = zeros(5);
       M(1,5) = g3 * c123;
       M(1,4) = g2*c12 + M(1,5);
       M(1,3) = g1*c1 + M(1,4);
       M(2,5) = g3 * s123;
       M(2,4) = g2 * s12 + M(2,5);
       M(2,3) = g1 * s1 + M(2,4);
       M(3,5) = k2 + m1*c23 + m2*c3;
       M(3,4) = k3 + m2*c3 + m3*c2 + M(3,5);
       M(4,5) = k2 + m2 * c3;
       % Construct the Mass Matrix:
       M = M + M' + diag(dM);
       
       dM = zeros(5,5,5);
       HM = zeros(5,5,5,5);
       % Calculate the first derivative of the mass matrix
       if nargout >= 2
           
           % Gradient wrt q5
           dM(1,3:5,5) = -g3 * s123;
           dM(2,3:5,5) = g3 * c123;
           dM(3,3:5,5) = -[(m1*s23 + m2*s3), m1 * s23 + 2*m2 * s3, m1 * s23 + m2 * s3];
           dM(4,4:5,5) = -[m2 * s3, m2 * s3];
           dM(:,:,5) = dM(:,:,5) + dM(:,:,5)';
           
           % Gradient wrt q4
           dM(1,:,4) = dM(1,:,5) - [0,0,g2*s12, g2*s12, 0];
           dM(2,:,4) = dM(2,:,5) + [0,0,g2*c12, g2*c12, 0];
           dM(3,:,4) = [0,0, -(m3*s2 + m1*s23), -(m1*s23 + m3*s2), -m1*s23];
           dM(:,:,4) = dM(:,:,4) + dM(:,:,4)';
           
           % Gradient wrt q3
           dM(1,:,3) = dM(1,:,4) - [0, 0, g1 * s1, 0, 0];
           dM(2,:,3) = dM(2,:,4) + [0, 0, g1 * c1, 0, 0];
           dM(:,:,3) = dM(:,:,3) + dM(:,:,3)';
           
           % Calculate the Hessian of the Mass Matrix
           if nargout >= 3
               
               % Only the Hessians related to the joint angles matter
               % Hessians with q5
               HM(1,:,5,3:5) = repmat([0, 0, -g3 * c123 * [1, 1, 1]], 1, 1, 1, 3);
               HM(2,:,5,3:5) = repmat([0, 0, -g3 * s123 * [1, 1, 1]], 1, 1, 1, 3);
               HM(3,:,5,5) = [0, 0, -(m1 * c23 + m2 * c3), -(m1 * c23 + 2*m2 * c3), -(m1*c23 + m2*c3)];
               HM(3,:,5,4) = [0, 0, -m1*c23, -m1*c23, -m1*c23];
               HM(4,:,5,5) = [0,0,0, -m2*c3, -m2*c3];
               
               % Hessians q4
               HM(:,:,4,5) = HM(:,:,5,4);
               HM(1:2,5,4,4) = HM(1:2,5,4,5);
               HM(1,3:4,4,4) = HM(1,3:4,5,4) - g2*c12;
               HM(2,3:4,4,4) = HM(2,3:4,5,4) - g2*s12;
               HM(3,3:5,4,4) = -[m3*c2 + m1*c23, m1*c23 + m3*c2, m1*c23];
               HM(1,3:5,4,3) = HM(1,3:5,5,3) - g2*c12*[1,1,0];
               HM(2,3:5,4,3) = HM(2,3:5,5,3) - g2*s12*[1,1,0];
                           
               % Hessians q3
               HM(:,:,3,5) = HM(:,:,5,3);
               HM(:,:,3,4) = HM(:,:,4,3);
               HM(1:2,:,3,3) = HM(1:2,:,4,3);
               HM(1,3,3,3) = HM(1,3,3,3) - g1*c1;
               HM(2,3,3,3) = HM(2,3,3,3) - g1*s1;
               
               % Symmetrize the Hessians
               HM = HM + permute(HM, [2,1,3,4]);
           end
       end
       
       
       
    end
    function [C, dC] = coriolisMatrix(self, q, dq)
            %% CORIOLISMATRIX: Matrix of Coriolis and Centripetal Effects
            %
            %   [C,dC] = coriolisMatrix(plant, q, dq) returns the matrix of
            %   Coriolis and Centripetal Effects, C, and its gradient, dC,
            %   for the PLANT model at configuration q and configuration
            %   rate dq. The Coriolis Matrix is used in the calculation of
            %   the PLANT dynamics:
            %
            %       M*ddq + C*dq + N = B*u
            %
            %   ARGUMENTS:
            %       PLANT:  a DIFFERENTIABLELAGRANGIAN Model
            %       q:      nQ x 1 vector of configuration variables
            %       dq:     nQ x 1 vector of configuration rates
            %
            %   RETURN VALUES:
            %       C:      nQ x nQ double matrix of Coriolis and
            %           centripetal effects
            %       dC:     nQ x nQ x 2nQ array of gradients for the 
            %           Coriolis matrix C. The first nQ pages are the
            %           gradients with respect to q, the second nQ pages
            %           are the gradients with respect to dq
            
            % Get the mass matrix gradient and Hessian to calculate the Coriolis Matrix
            [M, dM, d2M] = self.massMatrix(q);
            % Calculate the chirstoffel symbols
            G = 0.5 * (dM + permute(dM, [1,3,2]) - permute(dM, [3,2,1]));
            % Calculate the partials of the Christoffel symbols
            dG = 0.5 * (d2M + permute(d2M,[1,3,2,4]) - permute(d2M, [3,2,1,4]));
            % Calculate the Coriolis Matrix and its gradient using tensor multiplication
            C = times(G,reshape(dq,1,1,self.numQ));
            C = squeeze(sum(C,3));
            dG = times(dG, reshape(dq,1,1,self.numQ,1));
            dC = squeeze(sum(dG,3));
            % The gradient wrt dq is the Christoffel Symbols, G
            dC = cat(3, dC, G);   
    end
    function [N, dN] = gravityMatrix(self, q)

        
        % Inertial Constants
        k1 = self.baseMass + sum(self.masses);
        g1 = (self.masses(1)*self.com(1) + self.masses(2) + self.masses(3))*self.lengths(1);
        g2 = (self.masses(2)*self.com(2) + self.masses(3))*self.lengths(2);
        g3 = self.masses(3)*(self.com(3) - self.heelFraction)*self.lengths(3);
        
        % Sines and Cosines
        s1 = sin(q(3));
        c1 = cos(q(3));
        s12 = sin(q(3) + q(4));
        s123 = sin(q(3) + q(4) + q(5));
        c12 = cos(q(3) + q(4));
        c123 = cos(q(3) + q(4) + q(5));
       
        
        % Gravity Matrix
        N = zeros(5,1);
        N(2) = k1;
        N(5) = g3*s123;
        N(4) = g2*s12 + N(5);
        N(3) = g1*s1 + N(4);
        
        N = 9.81 * N;
        
        % Gradient of the gravitational effects
        if nargout == 2
           dN = zeros(5);
           
           dN(3:5,5) = g3*c123;
           dN(3:5,4) = dN(3:5,5) + g2*c12*[1;1;0];
           dN(3:5,3) = dN(3:5,4);
           dN(3,3) = dN(3,3) + g1*c1;
           
           dN = dN * 9.81;
        end
        
    end
    function [B, dB] = controllerMatrix(~,~)
      
       B = [zeros(2,3);eye(3)];
       if nargout == 2
           dB = zeros(5,3,5);
       end
    end
    %% -------------- VISUALIZATION METHODS -------------------- %%
    function [x,y] = draw(self, q, ax) 
        % First get the positions of all the joint centers
        [x,y] = self.positions(q);
        % Calculate the corners of the block
        th = -pi:0.01:pi;
        base_x = x(1) + self.base_radius * cos(th);
        base_y = y(1) + self.base_radius * sin(th);
        % Create a single data vector for plotting
        x = [base_x, nan, x];
        y = [base_y, nan, y];
        % Set x and y limits
        xrange = sum(self.lengths) + self.base_radius;
        yrange = sum(self.lengths) + self.base_radius;
        xlims = [x(1)-xrange, x(1)+xrange];
        ylims = [y(1)-yrange, y(1)+yrange];
        % Plot the data
        if nargout == 0
            % Create the floor
            [xfloor, yfloor] = self.terrain.draw(xlims, 100);
            x = [x,nan,xfloor];
            y = [y,nan,yfloor];
            if nargin == 3
                plot(ax,x,y,'LineWidth',2);
            else
                figure();
                plot(x,y,'LineWidth',2);
            end           
            xlim(xlims);
            ylim(ylims);
        end
    end
    function [x,y] = visualize(self, q, ax)
       %% VISUALIZE: Creates a visualization of the hopper at a specific configuration 
        
        
       if nargin == 2 && nargout == 0
           figure();
           ax = gca;
       end
       %% Do kinematics to get the joint positions
       [x,y] = self.positions(q);
       
       %% Get linkage drawings for each of the links
       % Base linkage (a ball)
       [x_base, y_base] = self.drawLinkage(x(1), y(1), 0, self.base_radius, 0);
       % Link 1 linkage
       [x_1, y_1] = self.drawLinkage(x(1), y(1), self.lengths(1), self.link_radius, 3*pi/2+ q(3));
       % Link 2 linkage
       [x_2, y_2] = self.drawLinkage(x(2), y(2), self.lengths(2), self.link_radius, 3*pi/2 + q(3) + q(4));
       % Link 3 linkage
       [x_3, y_3] = self.drawLinkage(x(5), y(5), self.lengths(3), self.link_radius, 3*pi/2 + q(3) + q(4) + q(5));
       
       %% Draw all the linkages
       if nargout == 0
           hold(ax, 'on');
           c = lines(4);
           patch(ax, x_1, y_1, c(2,:));
           patch(ax, x_2, y_2, c(3,:));
           patch(ax, x_3, y_3, c(4,:));
           % Draw the base linkage last.
           patch(ax, x_base, y_base, c(1,:));
           axis equal;
           yl = ylim;
           % Draw the terrain
           [xterrain, yterrain] = self.terrain.draw(xlim, 100);
           xterrain = [xterrain(:); xterrain(end); xterrain(1); xterrain(1)];
           yterrain = [yterrain(:); yl(1); yl(1); yterrain(1)];
           patch(ax, xterrain, yterrain, [0.4,0.4,0.4]);
       else
           x = [x_1(:), x_2(:), x_3(:), x_base(:)];
           y = [y_1(:), y_2(:), y_3(:), y_base(:)];
       end
    end
    function q = inverseKinematics(obj, xb, xf, phi)
        %% INVERSE KINEMATICS: Calculates configuration variables from cartesian positions
        %
        %    q = inverseKinematics(MODEL, X) calculates the
        %    configuration variables q for the system MODEl that achieve
        %    the desired end-effector cartesian position X.
        %
        %    ARGUMENTS:
        %        MODEL: A ContactDrivenCart object
        %        X: 4x1 double of Cartesian block and end-effector positions
        %        Q: 3x1 double of configuration coordinates
        
        % Calculate heel position
        xh = xf(1) - (1 - obj.heelFraction)*obj.lengths(3)*cos(phi);
        yh = xf(2) - (1 - obj.heelFraction)*obj.lengths(3)*sin(phi);
        % Calculate hip and knee angles
        % Relative heel position
        r = [xh - xb(1); yh - xb(2)];
        % Knee angle
        q2 = -acos((sum(r.^2) - sum(obj.lengths(1:2).^2))./(2*obj.lengths(1)*obj.lengths(2)));
        % Hip angle
        beta = atan2(r(1), r(2));
        alpha = atan2(obj.lengths(2)*sin(q2), obj.lengths(1) + obj.lengths(2)*cos(q2));
        q1 = pi - (beta + alpha);
        
        % Calculate ankle angle
        q3 = phi - q1 - q2 + pi/2;
        
        % Wrap the angles between -pi and pi
        q = wrapToPi([q1, q2, q3]);
        
        % Collect the configuration
        q = [xb(1), xb(2), q];
    end
end
end