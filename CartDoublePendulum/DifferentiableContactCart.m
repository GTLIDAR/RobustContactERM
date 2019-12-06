classdef DifferentiableContactCart < DifferentiableContactDynamics & DifferentiableLagrangian
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties
        % Parameters for the cart
        cartMass = 1;
        cartHeight = 1;
        % Parameters of the rods
        masses = [1,1];
        lengths = [1,1];
        com = [0.5, 0.5];
        inertias = [1, 1];
        % Visualization parameters
        cart_width = 1;
        cart_height = 1;
    end
    
    properties (Dependent)
        massFunctional;
        potentialFunctional;
        positionFunctional;
    end
    
    methods
        function obj = DifferentiableContactCart(terrain_geom)
            obj = obj@DifferentiableLagrangian(3,2);            
            if nargin == 1
               obj.terrain = terrain_geom; 
            end
        end
        
        function [w, dw, ddw] = configurationTransformation(~, q)
            
            % Transform the configuration into a set of coordinates that
            % includes all the information for calculating the dynamics and
            % their derivatives.
            w = [1, q(1), sin(q(2)), sin(q(2) + q(3)), cos(q(2)), cos(q(2)+q(3)), sin(q(3)), cos(q(3))]';
            % Partial derivative operators (sparse)
            d_dq1 = sparse(2,1,1,8,8);
            d_dq2 = sparse([3,4,5,6],[5,6,3,4],[1,1,-1,-1],8,8);
            d_dq3 = sparse([4,6,7,8],[6,4,8,7],[1,-1,1,-1],8,8);
            % Calculate the first and second mixed partial derivatives of 
            % the transformation
            dw = [d_dq1*w, d_dq2*w, d_dq3*w];
            ddw = [d_dq1 * dw, d_dq2*dw, d_dq3 * dw];
        end
        function [B, dB] = controllerMatrix(~, ~)
            B = [0,0;
                1,0;
                0,1];
            dB = zeros(3,2,3);
        end
        function fM = get.massFunctional(obj)
        
            % Calculate the individual constant terms in the mass matrix
            mTotal = sum(obj.masses) + obj.cartMass;
            g1 = obj.lengths(1) * (obj.masses(2) + obj.com(1)*obj.masses(1));
            g2 = obj.lengths(2) * obj.masses(2) * obj.com(2);            
            g3 = obj.masses(2) * (obj.com(2)*obj.lengths(2))^2;
            g4 = obj.masses(1) * (obj.com(1) * obj.lengths(1))^2 + obj.inertias(1) + obj.masses(2)*obj.lengths(1)^2 + g3;
            g5 = obj.masses(2) * obj.com(2) *obj.lengths(1) * obj.lengths(2);
            g6 = obj.inertias(2) + g3;
            
            
            % Collect the terms and their indices
            mvec = [mTotal, g1, g2, g2, g4, 2*g5, g6, g5, g6];
            rIdx = [1, 2, 2, 3, 4, 4, 5, 5, 6];
            cIdx = [1, 5, 6, 6, 1, 8, 1, 8, 1];
            % Set up the functional as a sparse matrix
            fM = sparse(rIdx, cIdx, mvec, 6,8);
        end

        function fV = get.potentialFunctional(obj)
            
            v_vec = 9.81 * [sum(obj.masses) + obj.cartMass, -(obj.com(1)*obj.masses(1) + obj.masses(2))*obj.lengths(1), ...
                - obj.com(2)*obj.masses(2)*obj.lengths(2)];
            fV = sparse([1, 1, 1], [1, 5, 6] ,v_vec, 1, 8);            
        end
        
        function fP = get.positionFunctional(obj)
        
            pVec = [1, obj.lengths, obj.cartHeight, -obj.lengths];
            fP = sparse([1,1,1,2,2,2], [2, 3, 4, 1, 5, 6], pVec, 2, 8);
        end
        
         %% -------------------- VISUALIZATION METHODS ----------------- %%
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
        
        function [x,y] = positions(obj,q)
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

