classdef PendulumDrivenCart < FreeCartDoublePendulum
  
    properties (Hidden)
        initialHeight(1,1) double = 1;
    end
    
methods
    function self = PendulumDrivenCart(varargin)
        self = self@FreeCartDoublePendulum(varargin{:});
    end
    %% --------------- KINEMATIC PROPERTIES -------------------  %%
    function [x,y] = positions(self,q)
       % Returns the positions of the links in the system
       q = [q(1);self.initialHeight;q(2:end)];
       [x,y] = positions@FreeCartDoublePendulum(self,q);
    end
    function J = jacobian(self, q)
       % Returns the endpoint Jacobian for the system
       % The endpoint Jacobian relates the generalized velocities to the
       % Cartesian coordinate velocities of the end of the second pendulum
       
       q = [q(1);0;q(2:end)];
       J = jacobian@FreeCartDoublePendulum(self,q);
       % Eliminate the ydot component
       J = J(:,[1,3:end]);
       
    end
    %% ---------------- DYNAMIC PROPERTIES --------------------- %%
    function M = inertiaMatrix(self, q)
        
        q = [q(1);0;q(2:end)];
        M = inertiaMatrix@FreeCartDoublePendulum(self,q);
        M = M([1,3:end],[1,3:end]);
    end
    function C = coriolisMatrix(self, q, qdot)
        q = [q(1);0;q(2:end)];
        qdot = [qdot(1);0;qdot(2:end)];
        C = coriolisMatrix@FreeCartDoublePendulum(self,q,qdot);
        C = C([1,3:end],[1,3:end]);
    end
    function N = gravityMatrix(self, q)
        % Returns gravity and conservative forces for the Pendulum Driven
        % Cart
        q = [q(1);0;q(2:end)];
        N = gravityMatrix@FreeCartDoublePendulum(self,q);
        N = N([1,3:end]);
    end
    function B = inputMatrix(~,~)
       B = [0,0;
            1,0;
            0,1];
    end
end
    
end