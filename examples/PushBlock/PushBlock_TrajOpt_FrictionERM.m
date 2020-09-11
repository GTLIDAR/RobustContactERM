function PushBlock_TrajOpt_FrictionERM()
%PUSHBLOCKTRAJOPTERM_CHANGINGFRICTION: Runs multiple trajectory
%optimization routines to evaluate the effect of changing the friction
%coefficient on the resulting trajectory
%
%

%   Luke Drnach
%   January 21, 2020

% Friction coefficient bounds
mu_low = 0.3;   % Minimum friction coefficient
mu_high = 0.7;  % Maximum friction coefficient
N_mu = 10;      % Number of total trials to run

% Generate a vector of friction coefficients
mu = linspace(mu_low, mu_high, N_mu);

% Create a base model for the plant
dt = 0.01;
plant = Block();
plant.timestep = dt;
% Add the FrictionERM solver
sig = 0.01;     % Friction uncertainty
plant.contactSolver =  UncertainFrictionERM(plant,0.5,sig);

% Specify the initial and final conditions
x0 = [0, plant.height/2, 0, 0]';
xf = [5, plant.height/2, 0, 0]';

% Generate the time vector
Tf = 1; %Final time
t_init = 0:plant.timestep:Tf;
N = numel(t_init);

% Generate the initial trajectory
x_init = zeros(length(x0),length(t_init));
for n = 1:length(x0)
    % Linear initial trajectory
   x_init(n,:) = linspace(x0(n),xf(n),N); 
end

% Create a structure for all the trajectories
trajData = struct('friction',[],'plant',[],'xtraj',[],'utraj',[],'z',[],'F',[],'info',[],'forces',[]);

% Force selector matrix
S = [1 0 0;
    0 1 -1];

% Run trajectory optimization with different mu and store the results
for n = 1:N_mu
    % Set the friction coefficient
    plant.terrain.friction_coeff = mu(n);
    plant.contactSolver.mu = mu(n);
    % Solve the trajectory optimization problem
    [xtraj, utraj, z, F, info] = optimizePushBlock(plant, t_init, x_init);
    % Store the data 
    trajData(n).friction = mu(n);
    trajData(n).plant = plant;
    trajData(n).xtraj = xtraj;
    trajData(n).utraj = utraj;
    trajData(n).z = z;
    trajData(n).F = F;
    trajData(n).info = info;
    % Calculate the ground reaction forces
    f = zeros(2,length(t_init));
    [~,x] = getPointsFromTrajectory(xtraj);
    [~,u] = getPointsFromTrajectory(utraj);
    for k = 1:length(t_init)
       f(:,k) = S*plant.contactForce(x(1:2,k),x(3:4,k),u(k)); 
    end
    trajData(n).forces = f;
    
end

% Save the data structure
save('PushBlock_MultiFrictionTraj_FrictionERM.mat','trajData');

% Plot the data all at once
plotAllData(trajData);

end

function [xtraj, utraj,z,F,info] = optimizePushBlock(plant, t_init, x_init)

% Get the final time and number of knot points
Tf = t_init(end);
N = length(t_init);
% Create a trajectory optimization problem
prob = DirtranTrajectoryOptimization(plant, N, Tf);
% Add the running cost
prob = prob.addRunningCost(@cost);
% Add the initial and final value constraints
prob = prob.addStateConstraint(ConstantConstraint(x_init(:,1)), 1);
prob = prob.addStateConstraint(ConstantConstraint(x_init(:,end)), N);
% Set options for the solver
prob = prob.setSolver('snopt');
prob = prob.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
prob = prob.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);
prob = prob.setSolverOptions('snopt','ScaleOption',1);
prob = prob.setSolverOptions('snopt','IterationsLimit',20000);
% Wrap the initial trajectory guess
traj_init.x = PPTrajectory(foh(t_init, x_init));
% Solve the problem
disp('Running Trajectory Optimization');
pause(0.2);
tic;
[xtraj, utraj, z, F, info] = prob.solveTraj(t_init, traj_init);
toc

end

function plotAllData(trajData)

f1 = figure();  %Trajectories
f2 = figure();  %Contact forces

% loop over all the data
for n = 1:length(trajData)
    % Get the points from the trajectory
    [t,x] = getPointsFromTrajectory(trajData(n).xtraj);
    [~,u] = getPointsFromTrajectory(trajData(n).utraj);
    lgdstr = sprintf('\\mu = %0.2f',trajData(n).friction);
    % Plot the state and control
    figure(f1);
    subplot(3,1,1);
    plot(t,x(1,:),'LineWidth',1.5,'DisplayName',lgdstr);
    hold on;
    ylabel('Position');
    subplot(3,1,2);
    plot(t,x(3,:),'LineWidth',1.5);
    hold on;
    ylabel('Velocity');
    subplot(3,1,3);
    plot(t,u,'LineWidth',1.5);
    hold on;
    ylabel('Control');
    xlabel('Time');
    % Plot the reaction forces
    figure(f2);
    subplot(2,1,1);
    plot(t,trajData(n).forces(1,:),'LineWidth',1.5,'DisplayName',lgdstr);
    hold on;
    ylabel('Normal');
    title('Contact Forces');
    subplot(2,1,2);
    plot(t,trajData(n).forces(2,:),'LineWidth',1.5);
    hold on;
    ylabel('Tangential');
    xlabel('Time');
end
figure(f1);
subplot(3,1,1);
legend show;
legend boxoff;
figure(f2);
subplot(2,1,1);
legend show;
legend boxoff;

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