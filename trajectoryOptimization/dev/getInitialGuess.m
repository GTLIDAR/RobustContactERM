function [x,u,lambda] = getInitialGuess(plant, x0, xf, N, method)
%GETINITIALGUESS Summary of this function goes here
%   Detailed explanation goes here


switch method
    case 'zero'
        [x,u,lambda] = guessZero(plant,N);
    case 'linear'
        [x,u,lambda] = guessLinear(plant,x0,xf,N);
    case 'static'
        [x,u,lambda] = guessStatic(plant,x0,xf,N);
end


end

function [x,u,lambda] = guessZero(plant,N)

x = zeros(plant.getNumStates(),N);
u = zeros(plant.getNumInputs(),N);
[Jn,Jt] = plant.contactJacobian(zeros(plant.getNumPositions(),1));

lambda = zeros(size(Jn,1)+size(Jt,1),N);

end
function [x,u,lambda] = guessLinear(plant, x0, xf, N)

[x,u,lambda] = guessZero(plant,N);
for n = 1:length(x0)
   x(n,:) = linspace(x0(n),xf(n),N); 
end

end
function [x,u,lambda] = guessStatic(plant,x0,xf,N)

[x,u,lambda] = guessLinear(plant,x0,xf,N);
B = plant.controllerMatrix(x(1:plant.getNumPositions(), 1));
active = logical(sum(B,2));
% Determine the forces and controls necessary to oppose gravity at each
% time point (assuming contact with the terrain)
for n = 1:N
    [Jn,Jt] = plant.contactJacobian(x(1:plant.getNumPositions(),n));
    J = [Jn; Jt];
    G = plant.gravityMatrix(x(1:plant.getNumPositions(),n));
    % Solve for the forces from the inactive degrees of freedom
    Jp = J(:,~active);
    lambda(:,n) = pinv(Jp')*G(~active);
    % Make sure all the contact forces are positive
    lambda(lambda(:,n)<0,n) =  0;
    % Solve for a holding controller
    u(:,n) = G(active) - J(:,active)'*lambda(:,n);
end


end