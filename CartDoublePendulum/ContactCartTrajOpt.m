function  ContactCartTrajOpt()
%CONTACTCARTTRAJOPT Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   December 2, 2019

% Add in the necessary paths
here = pwd;
cd ..
add_drake;
cd(here);
addpath(genpath('PATH_LCP'));

name = 'ContactCart_LCP_TrajOpt';
% Create the plant model
%terrain = FlatTerrain();
%plant = DifferentiableContactCart(terrain);
plant = ContactDrivenCart();
plant.cartHeight = 1.5;
% Specify the initial and final conditions
x0 = [0, pi/3, -2*pi/3, 0, 0, 0]';    
xf = [10, pi/3, -2*pi/3, 0, 0, 0]'; 

% Specify the number of collocation points
N = 101;

% Specify the final time of the problem
Tf = 5;     % Five seconds for this problem

plant.timestep = Tf/(N-1);

% Create a Trajectory Optimization Problem 
% Note that contact is implicit in the dynamics
prob = DirtranTrajectoryOptimization(plant, N, Tf);

% Add a running cost
prob = prob.addRunningCost(@cost);
% Add a terminal cost
prob = prob.addFinalCost(@terminalCost);
% Add the initial and final value constraints
prob = prob.addStateConstraint(ConstantConstraint(x0),1);
prob = prob.addStateConstraint(ConstantConstraint(xf),N);

% Set the options for the solver
prob = prob.setSolver('snopt');
prob = prob.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-5);
prob = prob.setSolverOptions('snopt','MajorOptimalityTolerance',1e-6);

% Create the initial guess at the solution
t_init = linspace(0, Tf, N);
x_init = zeros(numel(x0), N);
for n = 1:numel(x0)
   x_init(n,:) = linspace(x0(n),xf(n),N); 
end
traj_init.x = PPTrajectory(foh(t_init,x_init));
% Solve the problem
disp("Running Trajectory Optimization");
tic;
[xtraj, utraj, z, F, info] = prob.solveTraj(t_init, traj_init);
toc

save(name,'xtraj','utraj');

% Notes on the output of prob.solveTraj:
%   Return values:
%       xtraj: The state trajectory as a piecewise polynomial
%       utraj: the control trajectory as a piecewise polynomial
%       z:     the decision variables for the problem
%       F:     the value of the objective function after running the solver
%       INFO:  the exit flag of the solve. INFO = 1 indicates a successful
%              solve
%       
% Visualize the results

% Get the results
[~,xt] = getPointsFromTrajectory(xtraj);
%[~,ut] = getPointsFromTrajectory(utraj);

% Get the configuration vector
q = xt(1:3,:);
figure();
ax = gca;

% Animate the simulation
draw = @(ax,x) plant.draw(x,ax);
%utilities.animator(ax,draw,x(1:3,:));
utilities.animator(ax,draw,q,[name,'.avi']);

end


function [g, dg] = cost(dt, x, u)
% Running cost
R  =eye(2);
g = 1/2 * u' * R * u;
% The differential cost
dg = [zeros(1,1+size(x,1)),u'*R];
end

function [h, dh] = terminalCost(t, x)
% Terminal Cost    
h = t;
% Differential Terminal Cost
dh = [1, zeros(1,size(x,1))];
end
