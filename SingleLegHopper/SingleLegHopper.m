classdef SingleLegHopper <  Manipulator & DifferentiableContactDynamics 
    
    properties
        baseMass (1,1) double {mustBePositive, mustBeFinite} = 10;
        masses(2,1) double {mustBePositive, mustBeFinite} = [1,1];
        lengths(2,1) double {mustBeNonnegative, mustBeFinite} = [1 1];
        com(2,1) double = [1/2, 1/2];
        inertias(2,1) double = [1/12, 1/12];
    end
    properties (Hidden)
        base_radius(1,1) double = 0.25;
        numQ = 4;
        numU = 2;
    end
    
methods
    function self = SingleLegHopper(varargin)
        %% SingleLegHopper: Class Constructor
        
        nQ = 4;
        nU = 2;
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
    function [M, dM, HM] = massMatrix(self, q)
       %% Returns the Generalized Mass Matrix for the System
       %
       %
       
       %% Calculate Inertial constants
       % Constants
       k1 = self.baseMass + sum(self.masses);
       k2 = self.inertias(2) + self.masses(2) * (self.com(2) * self.lengths(2))^2;
       k3 = self.inertias(1) + (self.masses(1)  * self.com(1)^2 + self.masses(2)) * self.lengths(1)^2;
       % Sine / Cosine Multipliers
       g1 = (self.masses(1)*self.com(1) + self.masses(2))*self.lengths(1);
       g2 = (self.masses(2)*self.com(2))*self.lengths(2);
       g3 = g2 * self.lengths(1);
       
       %% Sine and Cosine
       s1 = sin(q(3));
       s2 = sin(q(4));
       c1 = cos(q(3));
       c2 = cos(q(4));
       s12 = sin(q(3) + q(4));
       c12 = cos(q(3) + q(4));
       
       %% Construct the Mass Matrix
       % Diagonals
       dM = [k1, k1, k2 + k3 + 2*(g3 * c2),k3];
       % Upper Triangle
       M = zeros(4);
       M(1,4) = g2*c12;
       M(1,3) = g1*c1 + M(1,4);
       M(2,4) = g2 * s12;
       M(2,3) = g1 * s1 + M(2,4);
       M(3,4) = k2 + g3*c2;
       % Construct the Mass Matrix:
       M = M + M' + diag(dM);
       
       dM = zeros(4,4,4);
       HM = zeros(4,4,4);
       % Calculate the first derivative of the mass matrix
       if nargout >= 2         
           % Gradient wrt q4
           dM(1,:,4) = -[0,0,g2*s12, g2*s12];
           dM(2,:,4) =  [0,0,g2*c12, g2*c12];
           dM(3,:,4) = [0,0, -g3*s2,  -g3*s2];
           dM(:,:,4) = dM(:,:,4) + dM(:,:,4)';
           
           % Gradient wrt q5
           dM(1,:,3) = dM(1,:,4) - [0, 0, g1 * s1, 0];
           dM(2,:,3) = dM(2,:,4) + [0, 0, g1 * c1, 0];
           dM(:,:,3) = dM(:,:,3) + dM(:,:,3)';
           
           % Calculate the Hessian of the Mass Matrix
           if nargout >= 3
               
               % Only the Hessians related to the joint angles matter
               % Hessians q4
               HM(1,3:4,4,4) = -g2*c12*[1,1];
               HM(2,3:4,4,4) = -g2*s12*[1,1];
               HM(3,3:4,4,4) = -g3*c2*[1,1];
               HM(1,3:4,4,3) = -g2*c12*[1,1];
               HM(2,3:4,4,3) = -g2*s12*[1,1];
                           
               % Hessians q3
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
        g1 = (self.masses(1)*self.com(1) + self.masses(2))*self.lengths(1);
        g2 = self.masses(2)*self.com(2)*self.lengths(2);
        
        % Sines and Cosines
        s1 = sin(q(3));
        c1 = cos(q(3));
        s12 = sin(q(3) + q(4));
        c12 = cos(q(3) + q(4));
       
        
        % Gravity Matrix
        N = zeros(4,1);
        N(2) = k1;
        N(4) = g2*s12;
        N(3) = g1*s1 + N(4);
        
        N = 9.81 * N;
        
        % Gradient of the gravitational effects
        if nargout == 2
           dN = zeros(4);
          
           dN(3:4,4) = g2*c12*[1;1];
           dN(3:4,3) = dN(3:4,4);
           dN(3,3) = dN(3,3) + g1*c1;
           
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
        block_x = x(1) + self.base_radius * cos(th);
        block_y = y(1) + self.base_radius * sin(th);
        % Create a single data vector for plotting
        x = [block_x, nan, x];
        y = [block_y, nan, y];
        % Set x and y limits
        xrange = sum(self.lengths) + self.base_radius;
        yrange = sum(self.lengths) + self.base_radius;
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
       m = self.baseMass + sum(self.masses); 
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