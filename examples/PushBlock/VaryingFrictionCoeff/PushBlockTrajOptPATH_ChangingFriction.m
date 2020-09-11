function PushBlockTrajOptPATH_ChangingFriction()
%PUSHBLOCKTRAJOPTPATH_CHANGINGFRICTION: Runs multiple trajectory
%optimization routines to evaluate the effect of changing the friction
%coefficient on the resulting trajectory
%
%

%   Luke Drnach
%   January 21, 2020

%   Updated January 22, 2020
%       Now uses a nominal trajectory from a previous optimization as a
%       warm-start

name = 'PathLCP';

% Friction coefficient bounds
mu_low = 0.3;   % Minimum friction coefficient
mu_high = 0.7;  % Maximum friction coefficient
N_mu = 5;      % Number of total trials to run

% Generate a vector of friction coefficients
mu = linspace(mu_low, mu_high, N_mu);
% Load the nominal case
nominal = load('WarmStart/PushBlock_TrajOpt_Nominal.mat');
plant = nominal.plant;

% Get the warm start trajectories
traj_init.x = nominal.xtraj;
traj_init.u = nominal.utraj;
t_init = nominal.t_init;

% Specify the initial and final conditions
x0 = nominal.x0;
xf = nominal.xf;

% Create a structure for all the trajectories
trajData = struct('friction',[],'plant',[],'xtraj',[],'utraj',[],'z',[],'F',[],'info',[],'forces',[]);

% Force selector matrix
S = [1 0 0;
    0 1 -1];

% Run trajectory optimization with different mu and store the results
for n = 1:N_mu
    % Set the friction coefficient
    plant.terrain.friction_coeff = mu(n);
    % Solve the trajectory optimization problem
    [xtraj, utraj, z, F, info] = optimizePushBlock(plant, t_init, x0, xf, traj_init);
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
save(['PushBlock_MultiFrictionTraj_', name,'.mat'],'trajData');

% Plot the data all at once
plotAllData(trajData,name);

end

