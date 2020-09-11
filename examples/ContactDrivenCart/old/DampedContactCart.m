classdef DampedContactCart < ContactDrivenCart
    %
    %
    
    %   Luke Drnach
    %   February 24, 2020
    
    properties
        damping = 0.7;
    end
    
    methods
        function [C, dC] = coriolisMatrix(obj, q, dq)
           % Get the coriolis matrix from the undamped cart
           [C, dC] = coriolisMatrix@ContactDrivenCart(obj, q, dq);
           % Add the damping
           C(1,1) = C(1,1) + obj.damping;
        end
    end
end

