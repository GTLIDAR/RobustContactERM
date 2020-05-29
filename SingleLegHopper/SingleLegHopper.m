classdef SingleLegHopper <  Manipulator & DifferentiableContactDynamics 
    
    properties
        blockMass (1,1) double {mustBePositive, mustBeFinite} = 10;
        masses(2,1) double {mustBePositive, mustBeFinite} = [1,1];
        lengths(2,1) double {mustBeNonnegative, mustBeFinite} = [1 1];
        com(2,1) double = [1/2, 1/2];
        inertias(2,1) double = [1, 1];
    end
    properties (Hidden)
        block_radius(1,1) double = 0.25;
    end
    
methods
    function self = SingleLegHopper(varargin)
        %% SingleLegHopper: Class Constructor
        
        nQ = 4;
        nU = 2;
        self = self@Manipulator(nQ, nU);
        
        self = self.setInputLimits(-1000,1000);
        
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
       
       x = zeros(1,3);
       y = zeros(1,3);
       % The position of the block is the first two coordinates
       x(1) = q(1);
       y(1) = q(2);
       % The positions of the links are determined by their angles
       x(2) = x(1) + self.lengths(1) * sin(q(3));
       y(2) = y(1) - self.lengths(1) * cos(q(3));
       x(3) = x(2) + self.lengths(2) * sin(q(3) + q(4)); 
       y(3) = y(2) - self.lengths(2) * cos(q(3) + q(4)); 
    end
    function x = kinematics(self, q)
       % Returns the position of the endpoint
        x = zeros(2,1);
        x(1) = q(1) + self.lengths(1) * sin(q(3)) + self.lengths(2) * sin(q(3) + q(4));
        x(2) = q(2) - self.lengths(1) * cos(q(3)) - self.lengths(2) * cos(q(3) + q(4));
    end
    
    function [J, dJ] = jacobian(self, q)
       % Returns the endpoint Jacobian for the system
       % The endpoint Jacobian relates the generalized velocities to the
       % Cartesian coordinate velocities of the end of the second pendulum
       
       J = zeros(2,4);
       J(1:2,1:2) = eye(2);
       J(1,4) = self.lengths(2)*cos(q(3)+q(4));
       J(2,4) = self.lengths(2)*sin(q(3)+q(4));
       J(1,3) = self.lengths(1)*cos(q(3)) + J(1,4);
       J(2,3) = self.lengths(1)*sin(q(3)) + J(2,4);
       
       if nargout == 2
          % The derivative of the Jacobian wrt each of the configuration
          % variables
          dJ = zeros(2,4,4);
          dJ(1,3:4,4) = -J(2,4);
          dJ(2,3:4,4) =  J(1,4);
          dJ(1,3:4,3) = -J(2,3:4);
          dJ(2,3:4,3) =  J(1,3:4);
       end
    end
    %% ---------------- DYNAMIC PROPERTIES --------------------- %%
    function [M, dM] = massMatrix(self, q)
        
        % Calculate some geometric constants
        mT = sum(self.masses) + self.blockMass;
        g1 = self.lengths(1)*(self.masses(2) + self.com(1)*self.masses(1));
        g2 = self.com(2)*self.lengths(2)*self.masses(2);
        g3 = self.masses(2)*(self.com(2)*self.lengths(2))^2;
        g4 = self.com(2)*self.lengths(2)*self.lengths(1)*self.masses(2);
        g5 = (self.masses(1)*self.com(1)^2 + self.masses(2))*self.lengths(1)^2;

        % Calculate the diagonals of the mass matrix
        diagM = [mT;
              mT;
              g5 + self.inertias(1) + g3 + 2*g4*cos(q(4));
              self.inertias(2) + g3];
        % Calculate the upper triangle of the mass matrix
        M = zeros(4);
        
        M(1,4) = g2*cos(q(3) + q(4));
        M(2,4) = g2*sin(q(3) + q(4));
        
        M(1,3) = g1*cos(q(3)) + M(1,4);
        M(2,3) = g1*sin(q(3)) + M(2,4);
       
        M(3,4) = self.inertias(2) + g3 + g4*cos(q(4));
        
        % Use symmetry to fill in the remaining elements and fill in the
        % diagonals
        M = M + M' + diag(diagM);
        
        % The derivatives of the mass matrix wrt the configuration
        % variables
        if nargout == 2
            dM = zeros(4,4,4);
            dM(1,3:4,3) = -M(2,3:4);
            dM(2,3:4,3) =  M(1,3:4);
            dM(:,:,3) = dM(:,:,3) + dM(:,:,3)';
            
            dM(1,3:4,4) = -M(2,4);
            dM(2,3:4,4) =  M(1,4);
            dM(3,3:4,4) = [-1,-1]*g4*sin(q(4));
            dM(:,:,4) = dM(:,:,4) + dM(:,:,4)';
        end
    end
    function [C, dC] = coriolisMatrix(self, q, qdot)
        C = zeros(4);
        
        % Pre-calculate some geometric constants
        g1 = self.lengths(1)*(self.masses(2) + self.com(1)*self.masses(1));
        g2 = self.com(2)*self.lengths(2)*self.masses(2);
        g4 = self.com(2)*self.lengths(2)*self.lengths(1)*self.masses(2);
        
        % Fill in the nonzero components of the Coriolis Matrix
        C(1,4) = - g2 * sin(q(3) + q(4)) * (qdot(3) + qdot(4));
        C(2,4) =   g2 * cos(q(3) + q(4)) * (qdot(3) + qdot(4));
        C(1,3) = - g1 * qdot(3) * sin(q(3)) + C(1,4); 
        C(2,3) =   g1 * qdot(3) * cos(q(3)) + C(2,4);
        C(3,3) = - g4 * sin(q(4)) * qdot(4);
        C(3,4) = - g4 * sin(q(4)) * (qdot(3) + qdot(4));
        C(4,3) =   g4 * sin(q(4)) * qdot(3);
        
        if nargout == 2
           dC = zeros(4,4,8); 
           % Deriv wrt theta_1 (q(3))
           dC(1,:,3) = -C(2,:);
           dC(2,:,3) = C(1,:);
           % Deriv wrt theta_2 (q(4))
           dC(1,3:4,4) = -C(2,4);
           dC(2,3:4,4) = C(1,4);
           dC(3,3:4,4) = -g4*cos(q(4))*[qdot(4), qdot(3)+qdot(4)];
           dC(4,3,4) = g4*cos(q(4))*qdot(3);
           % Deriv wrt theta_2_dot (dq(4))
           dC(1,3:4,8) = -g2 * sin(q(3)+q(4));
           dC(2,3:4,8) = g2 * cos(q(3) + q(4));
           dC(3,3:4,8) = -g4*sin(q(4));
           % Deriv wrt theta_1_dot (dq(3))
           dC(1:3,4,7) = dC(1:3,4,8);
           dC(1:2,3,7) = dC(1:2,3,8) + g1 *[-sin(q(3)); cos(q(3))];
           dC(4,3,7) = -dC(3,3,8);
        end
    end
    function [N, dN] = gravityMatrix(self, q)
        % Returns gravity and conservative forces for the Pendulum Driven
        % Cart
        g1 = self.lengths(1)*(self.masses(2) + self.com(1)*self.masses(1));
        g2 = self.com(2)*self.lengths(2)*self.masses(2);
        
        N = zeros(4,1);
        N(2) = sum(self.masses) + self.blockMass;
        N(4) = g2*sin(q(3)+q(4)); 
        N(3) = g1*sin(q(3)) + N(4);
        
        N = N * 9.81;
        
        if nargout == 2
           dN = zeros(4,4); 
           dN(3:4,3:4) = g2*cos(q(3)+q(4));
           dN(3,3) = dN(3,3) + g1*cos(q(3));
           dN = dN * 9.81;
        end
    end
    function [B, dB] = controllerMatrix(~,~)
      
       B = [zeros(2);eye(2)];
       if nargout == 2
           dB = zeros(4,2,4);
       end
    end
    %% -------------- VISUALIZATION METHODS -------------------- %%
    function [x,y] = draw(self, q, ax) 
        % First get the positions of all the joint centers
        [x,y] = self.positions(q);
        % Calculate the corners of the block
        th = -pi:0.01:pi;
        block_x = x(1) + self.block_radius * cos(th);
        block_y = y(1) + self.block_radius * sin(th);
        % Create a single data vector for plotting
        x = [block_x, nan, x];
        y = [block_y, nan, y];
        % Set x and y limits
        xrange = sum(self.lengths) + self.block_radius;
        yrange = sum(self.lengths) + self.block_radius;
        xlims = [x(1)-xrange,x(1)+xrange];
        ylims = [y(1)-yrange, y(1) + yrange];
        % Create the floor
        x = [x];
        y = [y];
        % Plot the data
        if nargout == 0
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
    function m = totalMass(self)
       m = self.blockMass + sum(self.masses); 
    end
    function q = inverseKinematics(obj, x)
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
        
        % First, move the block to the desired X location
        q = zeros(4,1);
        q(1) = x(1);
        q(2) = x(2);
        % Now calculate the end effector position relative to the cart
        % position
        r = [x(3) - x(1);x(4) - x(2)];
        % Calculate the second joint angle
        % Note there are actually two solutions for the second joint
        % angle, but acos will only return one
        q(4) = -acos((sum(r.^2) - sum(obj.lengths.^2))/(2*obj.lengths(1)*obj.lengths(2)));
        % Calculate the first joint angle using atan2 to get a unique
        % solution
        beta = atan2(r(1), r(2));
        alpha = atan2(obj.lengths(2)*sin(q(4)),  obj.lengths(1) + obj.lengths(2)*cos(q(4)));
        q(3) = pi - (beta + alpha);
    end
end
end