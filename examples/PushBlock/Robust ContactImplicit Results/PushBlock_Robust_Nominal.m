function PushBlock_Robust_Nominal()
%PushBlock_TrajOpt Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   December 17, 2019

% Add in the necessary paths
here = pwd;
cd ..
add_drake;
cd(here);
addpath(genpath('PATH_LCP'));

name = 'PushBlock_Robust_Nominal';
% Create the plant model
dt = 0.01;
plant = Block();
plant.timestep = dt;
plant.terrain.friction_coeff = 0.5;
% Specify the initial and final conditions
x0 = [0,plant.height/2,0,0]';     
xf = [5,plant.height/2,0,0]';    

% Specify the final time of the problem
Tf = 1;     % One second for this problem
    
t_init = 0:plant.timestep:Tf;
N = numel(t_init);

% Create a Trajectory Optimization Problem 
% Note that contact is implicit in the dynamics
options.integration_method = RobustContactImplicitTrajectoryOptimizer.BACKWARD_EULER;
options.contactCostMultiplier = 5;
options.uncertainty_source = RobustContactImplicitTrajectoryOptimizer.NO_UNCERTAINTY;
options.complementarity_solver = RobustContactImplicitTrajectoryOptimizer.NONE;
options.frictionVariance = 0.1;
prob = RobustContactImplicitTrajectoryOptimizer(plant, N, Tf, options);

% Add a running cost
%prob = prob.addRunningCost(@(t,x,u)cost2(t,x,u,xf));
prob = prob.addRunningCost(@cost);
% Add a terminal cost
%prob = prob.addFinalCost(@terminalCost);
% Add the initial and final value constraints
prob = prob.addStateConstraint(ConstantConstraint(x0),1);
prob = prob.addStateConstraint(ConstantConstraint(xf),N);

% Set the options for the solver
prob = prob.setSolver('snopt');
prob = prob.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-6);     %Default: 1e-6
prob = prob.setSolverOptions('snopt','MajorOptimalityTolerance',1e-6);      %Default: 1e-6
prob = prob.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-3);     %Default: 1e-6
prob = prob.setSolverOptions('snopt','ScaleOption',1);                      %Default: 1
prob = prob.setSolverOptions('snopt','IterationsLimit',50000);              %Default: 10,000
prob = prob.setSolverOptions('snopt','ElasticWeight',10^4);                 %Default: 10^4
prob = prob.setSolverOptions('snopt','MajorIterationsLimit',5000);          %Default: 1000
prob = prob.setSolverOptions('snopt','SuperbasicsLimit',1500);
% Create the initial guess at the solution
%t_init = linspace(0, Tf, N);
x_init = zeros(4,length(t_init));
for n = 1:4
    x_init(n,:) = linspace(x0(n),xf(n),N);
end
traj_init.x = PPTrajectory(foh(t_init,x_init));
% Solve the problem
disp("Running Trajectory Optimization");
pause(0.5);
tic;
[xtraj, utraj, ltraj, z, F, info] = prob.solveTraj(t_init, traj_init);
toc

save(name,'plant','xtraj','utraj','ltraj','t_init','z','F','info','plant','x0','xf');

% Notes on the output of prob.solveTraj:
%   Return values:
%       xtraj: The state trajectory as a piecewise polynomial
%       utraj: the control trajectory as a piecewise polynomial
%       z:     the decision variables for the problem
%       F:     the value of the objective function after running the solver
%       INFO:  the exit flag of the solve. INFO = 1 indicates a successful
%              solve
%       
%% Visualize the results

% Get the results
[t,x] = getPointsFromTrajectory(xtraj);
[~,u] = getPointsFromTrajectory(utraj);
[~,l] = getPointsFromTrajectory(ltraj);

% Get the contact forces from the lambda trajectory, ltraj
S = [1,0,0,0;0,1,-1,0];
f = S*l;

% Plot the results
figure();
subplot(3,1,1);
plot(t,x(1,:));
ylabel('Position');
subplot(3,1,2);
plot(t,x(3,:));
ylabel('Velocity');
subplot(3,1,3);
plot(t(1:end-1),u(1:end-1));
ylabel('Control');
xlabel('Time (s)');

% Contact impulses
figure();
subplot(2,1,1);
plot(t(1:end-1),f(1,1:end-1));
ylabel('Normal');
title('Contact Forces');
subplot(2,1,2);
plot(t(1:end-1),f(2,1:end-1));
ylabel('Tangential');
xlabel('Time (s)');


% Get the configuration vector
q = x(1:3,:);

%utilities.animator(ax,draw,x(1:3,:));
utilities.trajectoryAnimator(plant, q, [], [name,'.avi']);
end


function [g, dg] = cost(dt, x, u)
% Running cost
R = eye(numel(u)); % Control weights
Q = eye(numel(x)); % State weights

g = 1/2 * (u' * R * u + x'*Q*x);

% The differential cost
dg = [zeros(1),x'*Q,u'*R];
end
function [g, dg] = cost2(~, x, u, xf)
% Running cost
R = eye(numel(u)); % Control weights
Q = eye(numel(x)); % State weights

g = 1/2 * ( (x - xf)'*Q*(x - xf) + u' * R * u);
% The differential cost
dg = [zeros(1),(x - xf)'*Q,u'*R];
end
