classdef Terrain2D
    % Two-dimensional terrains
    methods (Abstract)
        % NEAREST: Given some point in 2D space xA, return the nearest
        % point X on the terrain. 
        x = nearest(obj, xA);
        % BASIS: Given some point on the terrain X, return a set of vectors
        % N and T, where N is normal to the terrain at X and T is
        % tangential to the terrain at X
        [N, T] = basis(obj, X)
    end
end