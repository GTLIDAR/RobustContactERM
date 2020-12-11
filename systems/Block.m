classdef Block < Manipulator & DifferentiableContactDynamics
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
            nQ = 2;
            nU = 1;
            obj = obj@Manipulator(nQ, nU);
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
        function [x,y] = visualize(obj, q)
           [x,y] = obj.draw(q); 
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
    methods (Static)
        function plotTrajectories(t, x, u, f, name)
            narginchk(4,5);
            % Get contact forces from solution trajectory
            S = [1,0,0,0;0,1,-1,0];
            f = S*f;
            % Plot the results
            if nargin< 5 || isempty(name)
                name = "Block";
            end
            % Horizontal trajectory
            figure('Name',[name,' Horizontal Trajectory']);
            subplot(3,1,1);
            plot(t,x(1,:));
            xlim([0, t(end)]);
            ylabel('Position');
            title('Horizontal');
            subplot(3,1,2);
            plot(t,x(3,:));
            xlim([0, t(end)]);
            ylabel('Velocity');
            subplot(3,1,3);
            plot(t(1:end-1),u(1:end-1));
            xlim([0, t(end)]);
            ylabel('Control');
            xlabel('Time (s)');
            % Vertical trajectory
            figure('Name',[name,' Vertical Trajectory']);
            subplot(2,1,1);
            plot(t,x(2,:));
            xlim([0, t(end)]);
            ylabel('Position');
            title('Vertical');
            subplot(2,1,2);
            plot(t,x(4,:));
            xlim([0, t(end)]);
            ylabel('Velocity');
            xlabel('Time (s)');
            % Contact impulses
            figure('Name',[name,' Reaction Forces']);
            subplot(2,1,1);
            plot(t(1:end-1),f(1,1:end-1));
            xlim([0, t(end)]);
            ylabel('Normal');
            title('Contact Forces');
            subplot(2,1,2);
            plot(t(1:end-1),f(2,1:end-1));
            xlim([0, t(end)]);
            ylabel('Tangential');
            xlabel('Time (s)'); 
        end
    end 
end

