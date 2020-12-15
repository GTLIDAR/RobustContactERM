function PushBlock_ERM_NoWarmstart()
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
snoptOptions = setSolverOptions([], ...
    'IterationsLimit',50000,...
    'MajorIterationsLimit',5000,...
    'SuperbasicsLimit',1500,...
    'MajorFeasibilityTolerance',1e-6,...
    'MajorOptimalityTolerance',1e-6,...
    'MinorFeasibilityTolerance',1e-6,...
    'ScaleOption',1);
optimOptions = setOptimizerOptions([],'nPoints',101,...
    'integration_method',RobustContactImplicitTrajectoryOptimizer.BACKWARD_EULER,...
     'nlcc_mode',1,...
     'duration',[1,1],...
     'display',false,...
     'uncertainty_source',RobustContactImplicitTrajectoryOptimizer.FRICTION_UNCERTAINTY,...
     'frictionVariance',0.01,...
     'frictionCostMultiplier',1e4,...
     'compl_slack',0);
 
%% Create the plant model
dt = 0.01;
plant = Block();
plant.timestep = dt;
plant.terrain.friction_coeff = 0.5;
%% Specify the initial and final conditions
x0 = [0,plant.height/2,0,0]';     
xf = [5,plant.height/2,0,0]';    

optimOptions.stateConstraints(1) = setStateConstraint(x0, 1, 1:4);
optimOptions.stateConstraints(2) = setStateConstraint(xf, optimOptions.nPoints, 1:4);

%% Add the running costs
optimOptions.runningCost(1).cost = @controlCost;
optimOptions.runningCost(1).name = 'ControlCost';
optimOptions.runningCost(2).cost = @(t,x,u)stateCost(t,x,u,xf);
optimOptions.runningCost(2).name = 'StateCost';
optimOptions.display = false;
%% Generate initial guess trajectory
guess.t = linspace(0, optimOptions.duration(2), optimOptions.nPoints);
x_init = zeros(4,length(guess.t));
for n = 1:4
    x_init(n,:) = linspace(x0(n),xf(n),optimOptions.nPoints);
end
guess.traj.x = PPTrajectory(foh(guess.t,x_init));

%% Solve trajectory optimization
init_guess = guess;
% Generate a vector of friction coefficient variances
uncertainty = logspace(log10(0.0001), log10(1), 9);
tolerance = [1e-6,1e-8];
source = pwd;
% Tag
tag = [pwd, '/', name];
for n = 1:length(uncertainty)
    guess = init_guess;
    target1 = sprintf('Uncertainty_%3.0E',uncertainty(n));
    mkdir(target1);
    cd(target1);
    for m = 1:length(tolerance)
        % Change the folder
        target2 = sprintf('Tolerance_%3.0E', tolerance(m));
        mkdir(target2);
        cd(target2);
        % Run the optimization
        optimOptions = setOptimizerOptions(optimOptions, 'frictionVariance',uncertainty(n));
        snoptOptions = setSolverOptions(snoptOptions, 'MajorOptimalityTolerance',tolerance(m),'MajorFeasibilityTolerance',tolerance(m));
        soln = optimizePlant(plant, guess, optimOptions, snoptOptions);
        % Save the results and figures in another folder
        save([name, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');
        % Print a report
        printReport(plant, soln, snoptOptions, optimOptions, tag);
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
        % Change the initial guess
        guess.traj.x = soln.xtraj;
        guess.traj.u = soln.utraj;
        guess.traj.lambda = soln.ltraj;
        % Revert to parent directory
        cd ..
    end
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
