classdef StepTerrain < Terrain2D

    %   Luke Drnach
    %   February 19, 2020
    
    properties
       friction_coeff = 0.5;
       terrain_height = 0.0;
       step_location = 0.0;
       step_height = 0.0;
    end
    
    methods
        function obj = StepTerrain()
        end
        function x = nearest(obj, xA)
            % Returns the nearest point on the terrain
            
            x = zeros(size(xA));
            x(1,:) = xA(1,:);
            x(2,:) = obj.terrain_height;
            x(2,xA(1,:) > obj.step_location) = obj.step_height;
            
        end
        function [N,T] = basis(obj, xB)
            
           nPoints = size(xB,2); 
           N = [0;1];
           T = [1;0];
           
           N = repmat(N, 1, nPoints);
           T = repmat(T, 1, nPoints);
        end
        function [x, y] = draw(obj, xlim, N)
           % Draws the terrain within the limits
           
           x = linspace(xlim(1), xlim(2), N);
           y = obj.terrain_height * ones(size(x));
           y(x > obj.step_location) = obj.step_height;
            
           % Find the point closes to the step location, and add an extra
           % point to draw the step correctly
           [~, idx] = min(abs(x - obj.step_location));
           if x(idx) == obj.step_location
               x = [x(1:idx), obj.step_location, x(idx+1:end)];
               y = [y(1:idx), obj.step_height, y(idx+1:end)];
           elseif x(idx) > obj.step_location
               x = [x(1:idx-1), obj.step_location, obj.step_location, x(idx:end)];
               y = [y(1:idx-1), obj.terrain_height, obj.step_height, y(idx:end)];
           else
               x = [x(1:idx), obj.step_location, obj.step_location, x(idx+1:end)];
               y = [y(1:idx), obj.terrain_height, obj.step_height, y(idx+1:end)];
           end
           if nargout == 0
              figure();
              plot(x, y, 'k-','LineWidth',2);
           end
        end
    end
end
