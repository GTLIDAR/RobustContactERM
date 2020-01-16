function PushBlock_TrajOpt()
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

name = 'PushBlock_TrajOpt';
% Create the plant model
dt = 0.01;
plant = Block();
plant.timestep = dt;

% Specify the initial and final conditions
x0 = [0,plant.height/2,0,0]';     
xf = [5,plant.height/2,0,0];    

% Specify the final time of the problem
Tf = 1;     % Five seconds for this problem

t_init = 0:plant.timestep:Tf;
N = numel(t_init);

% Create a Trajectory Optimization Problem 
% Note that contact is implicit in the dynamics
prob = DirtranTrajectoryOptimization(plant, N, Tf);

% Add a running cost
prob = prob.addRunningCost(@cost);
% Add a terminal cost
%prob = prob.addFinalCost(@terminalCost);
% Add the initial and final value constraints
prob = prob.addStateConstraint(ConstantConstraint(x0),1);
prob = prob.addStateConstraint(ConstantConstraint(xf),N);

% Set the options for the solver
prob = prob.setSolver('snopt');
prob = prob.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-6);
prob = prob.setSolverOptions('snopt','MajorOptimalityTolerance',1e-6);
prob = prob.setSolverOptions('snopt','ScaleOption',1);
prob = prob.setSolverOptions('snopt','IterationsLimit',20000);
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
[xtraj, utraj, z, F, info] = prob.solveTraj(t_init, traj_init);
toc

save(name,'plant','xtraj','utraj','z','F','info');

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

% Recalculate the contact forces
f0 = plant.contactForce(x(1:2,1),x(3:4,1),u(:,1));
f = zeros(length(f0),length(t));
f(:,1) = f0;
for n = 2:length(t)
   f(:,n) = plant.contactForce(x(1:2,n),x(3:4,n),u(:,n));
end

% Get the configuration vector
q = x(1:3,:);

%utilities.animator(ax,draw,x(1:3,:));
utilities.trajectoryAnimator(plant, q, [], [name,'.avi']);

% Plot the results
figure();
subplot(3,1,1);
plot(t,x(1,:));
ylabel('Position');
subplot(3,1,2);
plot(t,x(3,:));
ylabel('Velocity');
subplot(3,1,3);
plot(t,u);
ylabel('Control');
xlabel('Time (s)');

% Contact impulses
figure();
subplot(2,1,1);
plot(t,f(1,:));
ylabel('Normal');
title('Contact Forces');
subplot(2,1,2);
plot(t,f(2,:) - f(3,:));
ylabel('Tangential');
xlabel('Time (s)');
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