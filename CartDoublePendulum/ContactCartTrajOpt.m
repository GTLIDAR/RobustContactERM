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

% Load the kinematic solution
x = load('ContactCart_KinematicSequence.mat');
x_init = x.x;

name = 'ContactCart_LCP_TrajOpt';
% Create the plant model
dt = 0.01;
plant = ContactDrivenCart();
plant.timestep = dt;
plant.cartHeight = 1.5;

% Specify the initial and final conditions
x0 = [0,0];     % initial condition in cartesian coordinates
xf = [10,0];    % final condition in cartesian coordinates
% Calculate initial and final states
q0 = plant.inverseKinematics(x0);
qF = plant.inverseKinematics(xf);

x0 = [q0;zeros(3,1)];
xf = [qF; zeros(3,1)];

% Specify the number of collocation points
N = 101;

% Specify the final time of the problem
Tf = 5;     % Five seconds for this problem

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
prob = prob.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
prob = prob.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);
prob = prob.setSolverOptions('snopt','ScaleOption',1);
prob = prob.setSolverOptions('snopt','IterationsLimit',20000);
% Create the initial guess at the solution
t_init = linspace(0, Tf, N);
traj_init.x = PPTrajectory(foh(t_init,x_init));
% Solve the problem
disp("Running Trajectory Optimization");
tic;
[xtraj, utraj, z, F, info] = prob.solveTraj(t_init, traj_init);
toc

save(name,'xtraj','utraj','z','F','info');

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
tspan = xtraj.getBreaks();
xt = zeros(size(x0,1),numel(tspan));
ut = zeros(2,numel(tspan));


for n = 1:numel(tspan)
   xt(:,n) = xtraj.eval(tspan(n));
   ut(:,n) = utraj.eval(tspan(n));
end

% Get the configuration vector
q = xt(1:3,:);

%utilities.animator(ax,draw,x(1:3,:));
utilities.trajectoryAnimator(plant, q, [], [name,'.avi']);

end


function [g, dg] = cost(dt, x, u)
% Running cost
R = eye(numel(u)); % Control weights
Q = eye(numel(x)); % State weights
Q(1,1) = 0;

g = 1/2 * (u' * R * u + x'*Q*x);

% The differential cost
dg = [zeros(1),x'*Q,u'*R];
end

function [h, dh] = terminalCost(t, x)
% Terminal Cost    
h = t;
% Differential Terminal Cost
dh = [1, zeros(1,size(x,1))];
end

