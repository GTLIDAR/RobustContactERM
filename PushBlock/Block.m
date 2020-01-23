classdef Block < DifferentiableContactDynamics
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mass = 1;
        width = 1;
        height = 1;
        numQ = 2;
        numU = 1;
    end
    
    methods
        function obj = Block()
            % When we subclass Manipulator, we must provide
            % the number of configuration variables nQ and the
            % number of controls nU
            obj.numQ = 2;
            obj.numU = 1;
        end
        function [M, dM] = massMatrix(obj,q)
            M = obj.mass * eye(2);
            dM = zeros(2,2,2);
        end
        function [C, dC] = coriolisMatrix(obj, q, dq)
           C = zeros(2);
           dC = zeros(2,2,4);
        end
        function [N, dN] = gravityMatrix(obj, q)
           N = [0; obj.mass * 9.81];
           dN = zeros(2,2);
        end
        function [J, dJ] = jacobian(obj, q)
           J = eye(2);
           dJ = zeros(2,2,2);
        end
        function [B,dB] = controllerMatrix(obj,q)
           B = [1;0];
           dB = zeros(2,1,2);
        end
        function x = kinematics(obj,q)
            
            x = zeros(2,1);
            x(1) = q(1);
            x(2) = q(2) - obj.height/2;
        end
        function [x,y] = draw(obj,q,ax)
           
            % Calculate the corners of the block
            x = q(1) + 1/2*obj.width*[-1, -1, 1, 1, -1];
            y = q(2) + 1/2*obj.height*[1, -1, -1, 1, 1];
            % Set x and y limits
            xlims = q(1) + 1/2*obj.width*[-1,1];
            ylims = q(2) + 1/2*obj.height*[-1,1];
            % Expand the limits by 10% and make them square
            xlims = xlims + xlims .* sign(xlims) .* [-0.05, 0.05];
            ylims = ylims + ylims .* sign(ylims) .* [-0.05, 0.05];
            xr = range(xlims);
            yr = range(ylims);
            if xr > yr
               ylims = ylims * xr/yr; 
            else
               xlims = xlims * yr/xr;
            end
            
            % Plot the data
            if nargout == 0
                % Create the floor
                [xfloor, yfloor] = obj.terrain.draw(xlims, 100);
                if nargin < 3
                    figure();
                    ax = gca;
                end
                plot(ax,xfloor, yfloor, 'k-','LineWidth',2);
                hold on;
                plot(ax,x,y,'LineWidth',2);
                hold off;
                xlim(xlims);
                ylim(ylims);
            end
        end
    end
    
end

