classdef ContactDrivenCart < Manipulator & DifferentiableContactDynamics
    %% CONTACTDRIVENCART: An example for using Contact-Implicit Trajectory Optimization
    %
    %   ContactDrivenCart implements a two-link pendulum attached to a cart
    %   free to move in the horizontal direction but not the vertical
    %   direction. The plant is actuated only at the joints of the
    %   pendulums, and not actuated at the cart (underactuated of degree
    %   1). ContactDrivenCart requires a terrain to move between two points
    %   on the horizontal axis. As such, ContactDrivenCart is a simple
    %   example for Contact-Implicit Trajectory Optimization, as it
    %   requires contact with the ground to meet the objectives.
    %
    
    %   To be compatible with Drake, ContactDrivenCart subclasses the
    %   Manipulator class, and implements the manipulatorDynamics method.
    
    %   For a ContactDrivenCart, the configuration q consists of the
    %   horizontal position of the cart, q(1), the angle of the first link
    %   relative to vertical down, q(2), and the angle of the second link
    %   relative to the first link, q(3). This pattern is used throughout
    %   the implementation, and an identical pattern is used for the joint
    %   rates, dq.
    
    %   Luke Drnach
    %   November 19, 2019
    
    
    properties
        % Parameters for the cart
        blockMass = 1;
        cartHeight = 1;
        % Parameters of the rods
        masses = [1,1];
        lengths = [1,1];
        inertias = [1,1];
        com = [0.5, 0.5];
        % Configuration and Control space dimensions
        numQ = 3;
        numU = 2;
        % Parameters for the visualization methods
        cart_height = 1;
        cart_width = 1;
    end
   %% ------ Constructor and Required Methods for Manipulator ------ %% 
    methods
        function obj = ContactDrivenCart()
            %% ContactDrivenCart: Class Constructor
            %
            %   ContactDrivenCart creates a new instance of the
            %   ContactDrivenCart class, and returns it in OBJ
            %
            %   Syntax:
            %       obj = ContactDrivenCart();
            %
            %   Arguments:
            %       None
            %   Return Values:
            %       OBJ: An instance of the ContactDrivenCart class
            %
            %   ContactDrivenCart is a subclass of the Drake Class
            %   Manipulator. When subclassing Manipulator, we must provide
            %   the number of configuration variables nQ and the number of
            %   controls nU to the constructor for Manipulator. To modify
            %   the existing code to work without Drake, simply take out
            %   the lines referring to the Manipulator Class.
            
            % Create a Manipulator
            nQ = 3;
            nU = 2;
            obj = obj@Manipulator(nQ, nU);
            % Optionally set control limits (uncomment to include control
            % limits)
%             obj = setInputLimits(obj,-15,15);
        end 

    end
   %% ------ Methods defining the dynamic properties  -------------- %%
    methods 
        function [M,dM] = massMatrix(obj, q)
            %% MassMatrix: Calculates the system Mass Matrix
            %
            %   massMatrix calculates the mass matrix of the system, and
            %   the derivative of the mass matrix wrt the configuration, at
            %   a given configuration q.
            %
            %   Syntax:
            %       [M, dM] = massMatrix(obj, q);
            %       [M, dM] = obj.massMatrix(q);
            %
            %   Arguments:
            %       OBJ: A ContactDrivenCart object
            %       q:   3x1 double, the configuration vector
            %       
            %   Return Values
            %       M:  3x3 double, the mass matrix
            %       dM: 3x3x3 double, the derivatives of the mass matrix
            %           wrt the configuration, q. Note that dM(:,:,k) is
            %           the derivative of M wrt q(k).
            
            % Initialize the mass matrix array
            M = zeros(3);
            % Calculate some mass-geometric constants
            mT = sum(obj.masses) + obj.blockMass;
            g1 = obj.lengths(1)*(obj.masses(2) + obj.com(1)*obj.masses(1));
            g2 = obj.com(2)*obj.lengths(2)*obj.masses(2);
            g3 = obj.masses(2)*(obj.com(2)*obj.lengths(2))^2;
            g4 = obj.com(2)*obj.lengths(2)*obj.lengths(1)*obj.masses(2);
            g5 = (obj.masses(1)*obj.com(1)^2 + obj.masses(2))*obj.lengths(1)^2;
            % Calculate the diagonals of the mass matrix
            dM = [mT;
                g5 + obj.inertias(1) + g3 + 2*g4*cos(q(3));
                obj.inertias(2) + g3];
            % Calculate the upper triangle of the mass matrix
            M(1,3) = g2*cos(q(2) + q(3));
            M(1,2) = g1*cos(q(2)) + M(1,3);            
            M(2,3) = obj.inertias(2) + g3 + g4*cos(q(3));
            % Use symmetry to fill in the remaining elements and fill in the
            % diagonals
            M = M + M' + diag(dM);
            
            % Gradient tensor of the Mass matrix
            if nargout == 2
               dM = zeros(3,3,3);
               % Gradient with respect to block position is 0
               % Gradient with respect to the first joint
               dM(1,3,2) = - g2 * sin(q(2) + q(3));
               dM(1,2,2) = - g1 * sin(q(2)) + dM(1,3,2);
               dM(:,:,2) = dM(:,:,2) + dM(:,:,2)';
               % Gradient with respect to the second joint
               dM(1,2:3,3) = dM(1,3,2);
               dM(2,2:3,3) = -g4 * sin(q(3));
               dM(:,:,3)   = dM(:,:,3) + dM(:,:,3)'; 
            end
                        
        end
        function [C, dC] = coriolisMatrix(obj, q, dq)
            %% CoriolisMatrix: Matrix of Coriolis and Centripetal effects
            %
            %   coriolisMatrix returns the matrix of Coriolis and
            %   centripetal effects for the system. The Coriolis effects
            %   are those effects in which products of dissimilar
            %   velocities appear, while the centripetal effects are those
            %   force effects for which a product of similar velocities (a
            %   velocity squared) appears.
            %
            %   The vector of Coriolis and centripetal force effects can be
            %   calculated by multiplying the matrix by the joint rates:
            %   C*dq
            %
            %   Syntax:
            %       [C, dC] = coriolisMatrix(obj, q, dq);
            %       [C, dC] = obj.coriolisMatrix(q, dq);
            %
            %   Arguments:
            %       Obj:    a ContactDrivenCart object
            %       q:      3x1 double, the configuration vector
            %       dq:     3x1 double, the configuration rate vector
            %
            %   Return Values:
            %       C:      3x3 double, the Coriolis matrix
            %       dC:     3x3x6 double, the derivatives of the Coriolis
            %               matrix wrt q and dq. The first three pages are
            %               the derivatives wrt q; the last three pages are
            %               the derivatives wrt dq.
            
            % Initialize the coriolis matrix
            C = zeros(3);
            % Pre-calculate some geometric constants
            c1 = obj.com(2)*obj.masses(2)*obj.lengths(2);
            c2 = c1 * obj.lengths(1);
            c3 = (obj.masses(2) + obj.com(1)*obj.masses(1))*obj.lengths(1);
            % Fill in the nonzero components of the Coriolis Matrix
            C(1,3) = - c1 * sin(q(2) + q(3)) * (dq(2) + dq(3));
            C(1,2) = - c3 * dq(2) * sin(q(2)) + C(1,3);
            C(2,2) = - c2 * sin(q(3)) * dq(3);
            C(2,3) = - c2 * sin(q(3)) * (dq(2) + dq(3));
            C(3,2) =   c2 * sin(q(3)) * dq(2);
            % Gradient tensor of the Coriolis effects
            if nargout == 2
                dC = zeros(3,3,6);
                % Note that the derivatives wrt to the cart position q(1)
                % and the cart velocity dq(1) are all 0.
                % Derivative wrt q(2), the first link angle
                dC(1,3,2) = - c1 * cos(q(2) + q(3)) *(dq(2)+dq(3));
                dC(1,2,2) = -c3 * cos(q(2))*dq(2) +dC(1,3,2);
                % Derivative wrt q(3), the second link angle
                dC(1,2:3,3) = dC(1,3,2);
                dC(2,2,3) = -c2 * cos(q(3)) * dq(3);
                dC(2,3,3) = - c2 * cos(q(3)) * (dq(2) + dq(3));
                dC(3,2,3) = c2 * cos(q(3)) * dq(2);
                % Derivative wrt dq(2), the first link angular velocity
                dC(1,3,5) = -c1 * sin(q(2) + q(3));
                dC(1,2,5) = -c3 * sin(q(2)) + dC(1,3,5);
                dC(2,3,5) = -c2 * sin(q(3));
                dC(3,2,5) = - dC(2,3,5);
                % Derivative wrt dq(3), the second link angular velocity
                dC(1,2:3,6) = -c1 * sin(q(2) + q(3));
                dC(2,2:3,6) = -c2 * sin(q(3));
            end
        end  
        function [N,dN] = gravityMatrix(obj,q)
            %% GravityMatrix: Vector of Conservative Force Effects
            %
            %   gravityMatrix returns the vector of conservative (mostly
            %   gravitational) force effects for the system, and optionally
            %   the derivatives of the force effects wrt the configuration
            %   q.
            %
            %   Syntax:
            %       [N, dN] = gravityMatrix(obj, q);
            %       [N, dN] = obj.gravityMatrix(q);
            %
            %   Arguments:
            %       OBJ:    A ContactDrivenCart object
            %       q:      3x1 double, the configuration vector
            %       
            %   Return Values:
            %       N:      3x1 double, the conservative force effects
            %       dN:     3x3 double, the derivatives of N wrt q

            % Initialize the force effects vector
            N = zeros(3,1);
            % Calculate some geometric constants
            g1 = obj.com(2)*obj.masses(2)*obj.lengths(2);
            g2 = (obj.com(1)*obj.masses(1) + obj.masses(2))*obj.lengths(1);
            % Fill in the values for the force effects
            N(3) = g1*sin(q(2)+q(3));
            N(2) = g2*sin(q(2)) + N(3);
            N = N * 9.81;
            
            % Gradient matrix of the conservative forces
            if nargout == 2
               dN = zeros(3, 3); 
               % Gradients
               dN(3,2) = g1*cos(q(2) + q(3));
               dN(2,2) = g2*cos(q(2)) + dN(3,2);
               dN(:,3) = dN(3,2) * [0;1;1];
               dN = dN * 9.81;
            end
        end
        function [B,dB] = controllerMatrix(obj,~)
            %%  ControllerMatrix: The control selection matrix
            %
            %   controllerMatrix returns the control selection matrix for
            %   the system. Optionally, controllerMatrix also returns the
            %   derivatives of the control selection matrix wrt the
            %   configuration.
            %
            %   Syntax:
            %       [B, dB] = controllerMatrix(obj, q);
            %       [B, dB] = obj.controllerMatrix(q);
            %
            %   Arguments:
            %       OBJ:    a ContactDrivenCart object
            %       q:      3x1 double, the configuration vector
            %
            %   Return Values:
            %       B:      3x2 double, the controller selection matrix
            %       dB:     3x2x3 double, the derivatives of the controller
            %               selection matrix wrt the configuration vector q
            
            % The controller selection matrix
            B = [0, 0;
                1, 0;
                0, 1];
            % The derivatives of the selection matrix
            if nargout == 2
                dB = zeros(3,2,3);
            end
        end
        function [J, dJ] = jacobian(obj, q)
            %%  Jacobian: The Jacobian matrix for the end of the double pendulum
            %
            %   JACOBIAN returns the Jacobian matrix (the matrix of partial
            %   derivatives) for the end of the second pendulum in the
            %   ContactDrivenCart, evaluated at the current configuration.
            %
            %   Syntax:
            %       [J, dJ] = jacobian(obj, q);
            %       [J, dJ] = obj.jacobian(q);
            %
            %   Arguments:
            %       OBJ:    a ContactDrivenCart object
            %       q:      3x1 double, the configuration vector
            %
            %   Return Values:
            %       J:      2x3 double, the Jacobian matrix
            %       dB:     2x3x3 double, the derivatives of the Jacobian
            %               matrix wrt the configuration q
            
            % Initialize the Jacobian
            J = zeros(2,3);
            % Fill in the nonzero values
            J(1,1) = 1;
            J(1,3) = obj.lengths(2)*cos(q(2)+q(3));
            J(2,3) = obj.lengths(2)*sin(q(2)+q(3));
            J(1,2) = obj.lengths(1)*cos(q(2)) + J(1,3);
            J(2,2) = obj.lengths(1)*sin(q(2)) + J(2,3);
            
            % Gradient tensor for the Jacobian
            if nargout == 2 
                dJ = zeros(2,3,3);
                dJ(:,:,3) = [0, -obj.lengths(2) * sin(q(2) + q(3)) * [1 ,1];
                             0,  obj.lengths(2) * cos(q(2) + q(3)) * [1,1]];
                dJ(:,:,2) = [0, -obj.lengths(1) * sin(q(2)), 0;
                             0,  obj.lengths(1) * cos(q(2)), 0];
                dJ(:,:,2) = dJ(:,:,2) + dJ(:,:,3);
            end
        end
    end
    %% ------   Other Methods ------ %%
    methods
        function CoM = centerOfMass(obj, q)
            %% CENTEROFMASS: Returns the center of mass of the system 
            %
            %   centerOfMass calculates the center of mass for the entire
            %   ContactDrivenCart and returns it in an [x,y] pair. The
            %   center of mass is calculated for a given configuration.
            %
            %   Syntax:
            %       CoM = centerOfMass(obj, q);
            %       CoM = obj.centerOfMass(q);
            %
            %   Arguments:
            %       Obj:    a ContactDrivenCart object
            %       q:      3x1 double, the configuration vector
            %
            %    Return values:
            %        CoM:   2x1 double, the center of mass coordinates. The
            %               first value is the horizontal coordinate, the 
            %               second is the vertical coordinate.
            
            % Calculate the positions of the centers of mass for each body
            % in the system
            % For the cart:
            cart = [q(1), obj.cartHeight];
            % For the first link:
            link1 = [q(1) + obj.lengths(1)*obj.com(1)*sin(q(2)), obj.cartHeight - obj.lengths(1)*obj.com(1)*cos(q(2))];
            % For the second link:
            link2 = [q(1) + obj.lengths(1)*sin(q(2)) + obj.lengths(2)*obj.com(2)*sin(q(2)+q(3)), ...
                obj.cartHeight - obj.lengths(1)*cos(q(2)) - obj.lengths(2)*obj.com(2)*sin(q(2)+q(3))];
            % Take the weighted average to get the center of mass
            CoM = cart * obj.blockMass + link1 * obj.masses(1) + link2 * obj.masses(2);
            CoM = CoM ./(obj.blockMass + obj.masses(1) + obj.masses(2));
        end
        
        function [E, T, V] = energy(obj, q, dq)
           %% ENERGY: The total mechanical energy of the system
           %
           %    ENERGY calculates the total energy, the kinetic energy, and
           %    the potential energy of the system at a given state
           %    (configuration and rate pair)
           %
           %    Syntax:
           %        [E, T, V] = energy(obj, q, dq);
           %        [E, T, V] = obj.energy(q, dq);
           %
           %   Arguments:
           %       Obj:    a ContactDrivenCart object
           %       q:      3x1 double, the configuration vector
           %       dq:     3x1 double, the configuration rate vector
           %
           %    Return values:
           %        E:      scalar double, the total mechanical energy
           %        T:      scalar double, the kinetic energy
           %        V:      scalar double, the potential energy
           
           % Calculate the kinetic energy
           M = obj.massMatrix(q);
           T = dq'*M*dq / 2;
           % Calculate the potential energy
           V_b = obj.blockMass * obj.cartHeight * 9.81;
           V1 = (obj.cartHeight - obj.com(1) * obj.lengths(1) * cos(q(2))) * obj.masses(1) * 9.81;
           V2 = (obj.cartHeight - obj.lengths(1) * cos(q(2)) - obj.lengths(2) * obj.com(2) * cos(q(2)+q(3))) * obj.masses(2) * 9.81;
           V = V_b + V1 + V2;
           % Total energy
           E = T + V;
        end
        function x = kinematics(obj, q)
            
           x = zeros(2,1);
           x(1) = q(1) + obj.lengths(1) * sin(q(2)) + obj.lengths(2) * sin(q(2)+q(3));
           x(2) = obj.cartHeight - obj.lengths(1) * cos(q(2)) - obj.lengths(2) * cos(q(2) + q(3));
        end
        function [x,y] = draw(obj, q, ax)
            % First get the positions of all the joint centers
            [x,y] = obj.positions(q);
            % Calculate the corners of the block
            block_x = x(1) + 1/2*obj.cart_width*[-1, -1, 1, 1, -1];
            block_y = y(1) + 1/2*obj.cart_height*[1, -1, -1, 1, 1];
            % Create a single data vector for plotting
            x = [block_x, nan, x];
            y = [block_y, nan, y];
            % Set x and y limits
            xrange = sum(obj.lengths) + obj.cart_width;
            yrange = sum(obj.lengths) + obj.cart_height;
            xlims = [x(1)-xrange,x(1)+xrange];
            ylims = [y(1)-yrange, y(1) + yrange];
            % Plot the data
            if nargout == 0
                % Create the floor 
                [xfloor, yfloor] = obj.terrain.draw(xlims, 100);
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
        function [x, y] = positions(obj,q)
            % Returns the positions of the links in the system
            x = zeros(1,3);
            y = zeros(1,3);
            % The position of the block is the first two coordinates
            x(1) = q(1);
            y(1) = obj.cartHeight;
            % The positions of the links are determined by their angles
            x(2) = x(1) + obj.lengths(1) * sin(q(2));
            y(2) = y(1) - obj.lengths(1) * cos(q(2));
            x(3) = x(2) + obj.lengths(2) * sin(q(2) + q(3));
            y(3) = y(2) - obj.lengths(2) * cos(q(2) + q(3));
        end
        function q = iterativeIK(obj, x0, xF)
            
            x0 = x0(:);
            xF = xF(:);
            h = 0.01;                   % Step Size
            % Get the initial configuration using normal IK
            q0 = obj.inverseKinematics(x0);
            % Calculate the cartesian displacement
            dX = xF - x0;
            q = zeros(numel(q0),round(norm(dX)/h));
            q(:,1) = q0;
            n = 1;
            while norm(dX) > h
                J = obj.jacobian(q(:,n));
                q(:,n+1) = q(:,n) + h * pinv(J) * dX;
                x = obj.kinematics(q(:,n+1));
                dX = xF - x;
                n = n+1;
            end
            q = q(:,1:n);
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
           %        X: 2x1 double of Cartesian end-effector positions
           %        Q: 3x1 double of configuration coordinates
           
           % First, move the block to the desired X location
           q = zeros(3,1);
           q(1) = x(1);
           % Now calculate the end effector position relative to the cart
           % position
           r = [0;x(2) - obj.cartHeight ];
           % Calculate the second joint angle
           % Note there are actually two solutions for the second joint
           % angle, but acos will only return one
           q(3) = -acos((sum(r.^2) - sum(obj.lengths.^2))/(2*obj.lengths(1)*obj.lengths(2)));
           % Calculate the first joint angle using atan2 to get a unique
           % solution
           beta = atan2(r(1), r(2));
           alpha = atan2(obj.lengths(2)*sin(q(3)),  obj.lengths(1) + obj.lengths(2)*cos(q(3)));
           q(2) = pi - (beta + alpha);
        end
    end
end

 