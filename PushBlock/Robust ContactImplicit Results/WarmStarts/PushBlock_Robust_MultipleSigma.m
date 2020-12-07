function PushBlock_Robust_MultipleSigma()
%PUSHBLOCK_ROBUST_MULTIPLIESIGMA Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   February 5, 2020


name = 'ERM_Gaussian_1e4';

% Sigma Bounds
sigma_low = 0.0001;   % Minimum friction coefficient
sigma_high = 1;       % Maximum friction coefficient
N_sigma = 9;         % Number of total trials to run

% Generate a vector of friction coefficient variances
sigma = logspace(log10(sigma_low), log10(sigma_high), N_sigma);

options.uncertainty_source = RobustContactImplicitTrajectoryOptimizer.FRICTION_UNCERTAINTY;
options.distribution = RobustContactImplicitTrajectoryOptimizer.GAUSSIAN;
%options.complementarity_solver = RobustContactImplicitTrajectoryOptimizer.NONE;
options.integration_method = RobustContactImplicitTrajectoryOptimizer.BACKWARD_EULER;
options.contactCostMultiplier = 10000;
% Load the nominal case
nominal = load([pwd,'/PushBlock_Robust_RefinedNominal.mat']);
plant = nominal.plant;
t_init = nominal.t_init;
plant.terrain.friction_coeff = 0.5;
% Specify the initial and final conditions
x0 = nominal.x0;
xf = nominal.xf;

% Get the warm start trajectories
traj_init.x = nominal.xtraj;
traj_init.u = nominal.utraj;
traj_init.lambda = nominal.ltraj;

% Create a structure for all the trajectories
trajData = struct('sigma',[],'plant',[],'xtraj',[],'utraj',[],'z',[],'F',[],'info',[],'ltraj',[]);

% Run trajectory optimization with different mu and store the results
for n = 1:length(sigma)
    % Set the friction coefficient variance
    options.frictionVariance = sigma(n);

    % Solve the trajectory optimization problem
    [xtraj, utraj, ltraj, z, F, info] = optimizeRobustPushBlock(plant, t_init, x0, xf, traj_init,options);
    % Store the data 
    trajData(n).sigma = sigma(n);
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
save(['PushBlock_MultiSigma_', name,'.mat'],'trajData');

% Plot the data all at once
plotAllData(trajData,name,'sigma',nominal);

end

