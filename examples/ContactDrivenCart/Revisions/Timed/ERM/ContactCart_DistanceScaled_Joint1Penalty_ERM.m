function ContactCart_DistanceScaled_Joint1Penalty_ERM()
%%
%
%

%   Luke Drnach
%   April 17, 2020

%% Load the previous results and make an identifier
% Load the warm start
here = pwd;
cd('../Optimal');
warmstart = load([here, '/../Optimal/ContactCart_DistanceScaled_Joint1Penalty_Optimal.mat']);
cd(here);
fprintf('Loaded warm start from %s\n', warmstart.tag);
% Create a name and tag for this run
stk = dbstack;
name = stk.name;
tag = [here, '/', name];

%% Set up the values of sigma to test
sigma_low = 1e-3;
sigma_high = 1;
sigma = logspace(log10(sigma_low), log10(sigma_high), 7);

%% Extract the plant and other necessary information from the previous result
plant = warmstart.plant;
% Get the solver and optimization options
optimOptions = warmstart.optimOptions;
snoptOptions = warmstart.snoptOptions;
% Get the warmstart / initial guess
guess.traj.x = warmstart.soln.xtraj;
guess.traj.u = warmstart.soln.utraj;
guess.traj.lambda = warmstart.soln.ltraj;
guess.t = warmstart.soln.t;

%% Add a running cost to the problem
xf = optimOptions.stateConstraints(2).value;
optimOptions.runningCost = @(t, x, u)running_cost(t, x, u, xf);

%% Change to an ERM problem
optimOptions = setOptimizerOptions(optimOptions, 'heightVariance',1e-3,...
    'uncertainty_source',RobustContactImplicitTrajectoryOptimizer.DISTANCE_UNCERTAINTY,...
    'contactCostMultiplier',1e5,'display',false);

%% Solve the problem
soln(length(sigma)) = warmstart.soln;
q = cell(1,length(sigma));
label = q;
for n = 1:length(sigma)
    fprintf('\nRun %d Sigma = %0.6e',n, sigma(n));
    optimOptions = setOptimizerOptions(optimOptions, 'heightVariance',sigma(n));
    soln(n)  = optimizePlant(plant, guess, optimOptions, snoptOptions);
    % Print a report
    filename = sprintf('run_%d_report.txt',n);
    printReport(plant, soln(n), snoptOptions, optimOptions, warmstart.tag,filename);
    % Store the configuration for animation
    x = soln(n).xtraj.eval(soln(n).t);
    q{n} = x(1:3,:);
    label{n} = sprintf('\\sigma = %0.3e', sigma(n));
end
for n = 1:length(sigma)
    soln(n).sigma = sigma(n);
end
%% Save the results and visualize
save([name, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');

%% Make a combined animation
animationUtilities.multiTrajectoryAnimator(plant, q, xf(1:3), [name,'.avi'], label)
end

function [g, dg] = running_cost(t, x, u, xf)
% Cost Weight matrices
R = diag([1,1]);
Q = diag([1, 100, 10, 1, 100, 10]);
% Penalize state deviations
dx = x - xf;
% The running cost function
g = 0.5 * (dx'*Q*dx + u'*R*u);
% The differential cost
dg = [0, dx'*Q, u'*R];
end