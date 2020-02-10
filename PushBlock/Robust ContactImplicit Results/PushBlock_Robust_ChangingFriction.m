function PushBlock_Robust_ChangingFriction(warmstart)
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

if nargin == 0
    warmstart = true;
end

name = 'Robust_ERM_Gaussian_ModerateSigma';

% Friction coefficient bounds
mu_low = 0.3;   % Minimum friction coefficient
mu_high = 0.7;  % Maximum friction coefficient
N_mu = 10;      % Number of total trials to run

% Generate a vector of friction coefficients
mu = linspace(mu_low, mu_high, N_mu);
options.uncertainty_source = RobustContactImplicitTrajectoryOptimizer.FRICTION_UNCERTAINTY;
options.distribution = RobustContactImplicitTrajectoryOptimizer.GAUSSIAN;
options.frictionVariance = 0.01;
options.integration_method = RobustContactImplicitTrajectoryOptimizer.BACKWARD_EULER; %Backward Euler Integration
options.contactCostMultiplier = 1000000;
options.complementarity_solver = RobustContactImplicitTrajectoryOptimizer.NONE;


% Load the plant and conditions from the warm start
nominal = load('WarmStarts/PushBlock_Robust_NCP_NominalRefined.mat');
plant = nominal.plant;
% Specify the initial and final conditions
x0 = nominal.x0;
xf = nominal.xf;
t_init = nominal.t_init;

if warmstart    
    % Get the warm start trajectories
    traj_init.x = nominal.xtraj;
    traj_init.u = nominal.utraj;
    traj_init.lambda = nominal.ltraj;
else
    % No warm start for the trajectory
    N = length(t_init);
    x = zeros(length(x0), N);
    for k = 1:length(x0)
       x(k,:) = linspace(x0(k), xf(k), N); 
    end
    traj_init.x = PPTrajectory(foh(t_init,x_init));
end


% Create a structure for all the trajectories
trajData = struct('friction',[],'plant',[],'xtraj',[],'utraj',[],'z',[],'F',[],'info',[],'ltraj',[]);

% Run trajectory optimization with different mu and store the results
for n = 1:length(mu)
    % Set the friction coefficient
    plant.terrain.friction_coeff = mu(n);
    % Solve the trajectory optimization problem
    [xtraj, utraj, ltraj, z, F, info] = optimizeRobustPushBlock(plant, t_init, x0, xf, traj_init,options);
    % Store the data 
    trajData(n).friction = mu(n);
    trajData(n).plant = plant;
    trajData(n).xtraj = xtraj;
    trajData(n).utraj = utraj;
    trajData(n).z = z;
    trajData(n).F = F;
    trajData(n).info = info;
    % Calculate the ground reaction forces
    trajData(n).ltraj = ltraj;
    
end

% Save the data structure
save(['PushBlock_MultiFrictionTraj_', name,'.mat'],'trajData');

% Plot the data all at once
plotAllData(trajData,name,'friction',nominal);

end

