function PushBlock_ERM_Warmstart()
%PushBlock_TrajOpt Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   December 3 2020

% Add in the necessary paths
% here = pwd;
% cd ..
% add_drake;
% cd(here);

stk = dbstack;
name = stk.name;
%% Load the warmstart from the Nominal solution
nominal = load([pwd, '/../Nominal/Tolerance_1E-08/PushBlock_Nominal.mat']);
disp(['Loaded from ',nominal.tag]);

snoptOptions = nominal.snoptOptions;
optimOptions = setOptimizerOptions(nominal.optimOptions,...
     'uncertainty_source',RobustContactImplicitTrajectoryOptimizer.FRICTION_UNCERTAINTY,...
     'frictionVariance',0.01,...
     'frictionCostMultiplier',1e4);
 
%% Get the plant and warmstart
plant = nominal.plant;
guess.traj.x = nominal.soln.xtraj;
guess.traj.u = nominal.soln.utraj;
guess.traj.lambda = nominal.soln.ltraj;
guess.t = nominal.soln.t;
 
% Final state
xf = optimOptions.stateConstraints(2).value;

%% Add the running costs
optimOptions.runningCost(1).cost = @controlCost;
optimOptions.runningCost(1).name = 'ControlCost';
optimOptions.runningCost(2).cost = @(t,x,u)stateCost(t,x,u,xf);
optimOptions.runningCost(2).name = 'StateCost';
optimOptions.display = false;
%% Solve trajectory optimization
% Generate a vector of friction coefficient variances
uncertainty = logspace(log10(0.0001), log10(1), 9);
source = pwd;
% Tag
tag = [pwd, '/', name];
for n = 1:length(uncertainty)
    % Change the folder
    target = ['Uncertainty_',sprintf('%3.0E',uncertainty(n))];
    mkdir(target);
    cd(target);
    % Run the optimization
    optimOptions = setOptimizerOptions(optimOptions, 'frictionVariance',uncertainty(n));
    soln = optimizePlant(plant, guess, optimOptions, snoptOptions);
    % Save the results and figures in another folder
    save([name, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');
    % Print a report
    printReport(plant, soln, snoptOptions, optimOptions, nominal.tag);
    % Visualize solution trajectories
    try
        visualizeSolutionTrajectories(soln, plant, name);
    catch
        disp('Animation failed');
    end
    % Visualize complementarity relationships
    visualizeComplementarity(soln, plant, name);
    % Save figures and close
    saveFigsAndClose();
    % Revert to parent directory
    cd(source);
end
end

function [g,dg] = controlCost(dt, x, u)
R = 100*eye(numel(u)); % Control weights
g = 1/2 * ( u' * R * u);
% The differential cost
dg = [zeros(1,numel(dt)+numel(x)),u'*R];
end
function [g, dg] = stateCost(dt, x, u, xf)
Q = eye(numel(x)); % State weights

g = 1/2 * ( (x - xf)'*Q*(x - xf));
% The differential cost
dg = [zeros(numel(1,dt)),(x - xf)'*Q,zeros(1,numel(u))];
end
