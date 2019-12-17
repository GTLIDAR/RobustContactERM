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

