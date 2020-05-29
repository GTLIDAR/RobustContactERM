classdef FlatTerrain < Terrain2D
    %FlatTerrain - 2D flat terrain
    %
    %   
    
    %   Luke Drnach
    %   November 20, 2019
    
    properties
        friction_coeff = 0.5;
        terrain_height = 0;
        basis_mult = 1;
    end
    
    methods
        function obj = FlatTerrain()
        end
        
        function x = nearest(obj, xA)
            % NEAREST: Returns the nearest point on the terrain
            % 
            %   For FLAT TERRAIN, the nearest point to (xA, yA) is (xA, 0)
            x = zeros(size(xA));
            x(1,:) = xA(1,:);
            x(2,:) = obj.terrain_height;
        end
        function [N, T] = basis(obj, xB)
           % BASIS: Returns the local coordinate basis for the terrain
           %
           %   For a flat terrain, the basis is simply the world
           %   coordinates
           %
           %    RETURN VALUE
           %        N: 2x1 double, the normal vector to the terrain
           %        T: 2x1 double, the tangent vector to the terrain
           %
           N = obj.basis_mult*[0;1];
           T = obj.basis_mult*[1;0];
           
           % Replicate for each contact point
           nPoints = size(xB, 2);
           N = repmat(N, 1, nPoints);
           T = repmat(T, 1, nPoints);
        end        
        function [x,y] = draw(obj, xlim, N)
           % DRAW: Draws the terrain within the limits
           %
           %
           
           x = linspace(xlim(1),xlim(2),N);
           y = obj.terrain_height * ones(size(x));
        end
    end
    
end

