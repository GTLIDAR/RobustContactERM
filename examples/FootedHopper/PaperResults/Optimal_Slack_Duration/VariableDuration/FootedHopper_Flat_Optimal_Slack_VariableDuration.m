function FootedHopper_Flat_Optimal_Slack_VariableDuration()
%OPTIMAL_LARGEFOOT Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   May 19, 2020

%% Load the previous results and make an identifier
% Load the warm start
here = pwd;
feasible = load([here, '/../../Feasible_Slack_Duration/VariableDuration/Slack_0E+00/FootedHopper_Flat_Feasible_Slack_VariableDuration.mat']);
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
guess.traj.slacks = feasible.soln.slacks;
guess.traj.jltraj = feasible.soln.jltraj;
%% Add a running cost to the problem
xf = optimOptions.stateConstraints(2).value;
%optimOptions.runningCost = @(t, x, u)running_cost(t, x, u, xf, plant);
optimOptions.runningCost(1).cost = @controlCost;
optimOptions.runningCost(1).name = 'ControlCost';
optimOptions.runningCost(2).cost = @(t,x,u)stateCost(t,x,u,xf);
optimOptions.runningCost(2).name = 'StateCost';

% Turn on cost displaying
optimOptions = setOptimizerOptions(optimOptions,'display',true);
%% Solve the problem
snoptOptions = setSolverOptions(snoptOptions, 'MajorIterationsLimit',1e4,'IterationsLimit',1e5,'ScaleOption',2);
soln  = optimizePlant(plant, guess, optimOptions, snoptOptions);

%% Save the results and visualize
save([name, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');

% Print a report
printReport(plant, soln, snoptOptions, optimOptions, feasible.tag);

% Visualizer
try
    visualizeFootedHopper(plant, soln, name);
catch
    disp('Animation failed');
end
savefig(figure(1),'CostIterations');

end

function [g,dg] = controlCost(t, x, u)

% Cost weight matrix
R = 0.01*diag([1,1,1]);
% Quadratic Control Cost
g = 0.5 * u'*R*u;
% Differential cost
dg = [zeros(1,1+numel(x)), u'*R];

end
function [g,dg] = stateCost(t,x,u,xf)

% Cost weight matrix
Q = diag([1, 10, 10, 100, 100, 1, 1, 1, 1, 1]);
% Penalize state deviations
dx = x - xf;
% Quadratic state cost
g = 0.5 * dx'*Q*dx;
% Differential Cost
dg = [0, dx'*Q, zeros(1,numel(u))];
end