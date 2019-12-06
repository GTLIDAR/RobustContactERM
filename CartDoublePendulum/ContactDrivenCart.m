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
        % Terrain
        %terrain = FlatTerrain();
        numQ = 3;
        numU = 2;
        % Visualization methods
        cart_height = 1;
        cart_width = 1;
    end
    
    methods
        function obj = ContactDrivenCart()
            
            % When we subclass Manipulator, we must provide
            % the number of configuration variables nQ and the 
            % number of controls nU
            nQ = 3;
            nU = 2;
            obj = obj@Manipulator(nQ, nU);
            % We can optionally set input limits for the acrobot
            obj = setInputLimits(obj,-10,10);
        end 
        function [H,C,B, dH, dC, dB] = manipulatorDynamics(obj,q,dq)
            
            [H, dH] = obj.massMatrix(q);
            [C, dC] = obj.coriolisMatrix(q,dq);
            [N, dN] = obj.gravityMatrix(q);
            [B, dB] = obj.controllerMatrix(q);
            
            % Gradient wrt q
            dC_q = times(dC(:,:,1:obj.numQ), reshape(dq,1,obj.numQ,1));
            dC_q = squeeze(sum(dC_q, 2));
            
            dC_q = dC_q + dN;
            % Gradient wrt dq
            dC_dq = times(dC(:,:,obj.numQ+1:end), reshape(dq, 1, obj.numQ, 1));
            dC_dq = squeeze(sum(dC_dq, 2));
            dC_dq = dC_dq + C;
            
            % Gradient of C*q + N
            dC = [dC_q, dC_dq]';
            
            % Combine the coriolis and gravitational effects
            C = C*dq + N;
            
            % Reshape dM and dB
            dH = reshape(dH, [size(dH, 1), size(dH,2) * size(dH, 3)])';
            dB = reshape(dB, [size(dB, 1), size(dB, 2) * size(dB, 3)])';
        end
    end
   
    methods 
        function [M,dM] = massMatrix(obj, q)
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
            M = zeros(3);
            
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
                
                dC(1,3,2) = - c1 * cos(q(2) + q(3)) *(dq(2)+dq(3));
                dC(1,2,2) = -c3 * cos(q(2))*dq(2) +dC(1,3,2);
                
                dC(1,2:3,3) = dC(1,3,2);
                dC(2,2,3) = -c2 * cos(q(3)) * dq(3);
                dC(2,3,3) = - c2 * cos(q(3)) * (dq(2) + dq(3));
                dC(3,2,3) = c2 * cos(q(3)) * dq(2);
                
                dC(1,3,5) = -c1 * sin(q(2) + q(3));
                dC(1,2,5) = -c3 * sin(q(2)) + dC(1,3,5);
                dC(2,3,5) = -c2 * sin(q(3));
                dC(3,2,5) = - dC(2,3,5);
                
                dC(1,2:3,6) = -c1 * sin(q(2) + q(3));
                dC(2,2:3,6) = -c2 * sin(q(3));
                
            end
        end
        
        function [N,dN] = gravityMatrix(obj,q)
            % Returns gravity and conservative forces for the Pendulum Driven
            % Cart
            g1 = obj.com(2)*obj.masses(2)*obj.lengths(2);
            g2 = (obj.com(1)*obj.masses(1) + obj.masses(2))*obj.lengths(1);
            
            N = zeros(3,1);
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
           B = [0, 0;
                1, 0;
                0, 1];
           if nargout == 2
              dB = zeros(3,2,3); 
           end
        end
        function [J, dJ] = jacobian(obj, q)
            % Returns the endpoint Jacobian for the system
            % The endpoint Jacobian relates the generalized velocities to the
            % Cartesian coordinate velocities of the end of the second pendulum
            
            J = zeros(2,3);
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
            % Create the floor -- THIS WORKS WITH ONLY FLAT TERRAIN
            xfloor= xlims;
            yfloor = [0,0];
            x = [x,nan,xfloor];
            y = [y,nan,yfloor];
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
    end
end

