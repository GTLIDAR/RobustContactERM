function ContactCart_DistanceScaled_Joint1Penalty_Optimal()
%%
%
%

%   Luke Drnach
%   February 28, 2020

%% Load the previous results and make an identifier
% Load the warm start
here = pwd;
feasible = load([here, '/../Feasible/ContactCart_DistanceScaled_Feasible.mat']);
fprintf('Loaded warm start from %s\n', feasible.tag);
% Create a name and tag for this run
stk = dbstack;
name = stk.name;
tag = [here, '/', name];

%% Extract the plant and other necessary information from the previous result
plant = feasible.plant;
% Get the solver and optimization options
optimOptions = feasible.optimOptions;
snoptOptions = feasible.snoptOptions;
% Get the warmstart / initial guess
guess.traj.x = feasible.soln.xtraj;
guess.traj.u = feasible.soln.utraj;
guess.traj.lambda = feasible.soln.ltraj;
guess.t = feasible.soln.t;

%% Add a running cost to the problem
xf = optimOptions.stateConstraints(2).value;
optimOptions.runningCost = @(t, x, u)running_cost(t, x, u, xf);
optimOptions = setOptimizerOptions(optimOptions,'display',false);
%% Solve the problem
soln  = optimizePlant(plant, guess, optimOptions, snoptOptions);

%% Save the results and visualize
save([name, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');

% Print a report
printReport(plant, soln, snoptOptions, optimOptions, feasible.tag);

% Visualizer
visualizeCartTrajectories(plant, soln, name);
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