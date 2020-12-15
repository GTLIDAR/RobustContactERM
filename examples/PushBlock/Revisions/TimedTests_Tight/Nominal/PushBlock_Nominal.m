function PushBlock_Nominal()
%PushBlock_TrajOpt Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   December 11 2020

% Add in the necessary paths
stk = dbstack;
name = stk.name;
%% Set the options for the solver
% Set low-level options for SNOPT
% NOTE: DRAKE does not allow us to set the low-level option ElasticWeight
% for SNOPT.
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
     'uncertainty_source',RobustContactImplicitTrajectoryOptimizer.NO_UNCERTAINTY,...
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
%% Generate initial guess trajectory
guess.t = linspace(0, optimOptions.duration(2), optimOptions.nPoints);
x_init = zeros(4,length(guess.t));
for n = 1:4
    x_init(n,:) = linspace(x0(n),xf(n),optimOptions.nPoints);
end
guess.traj.x = PPTrajectory(foh(guess.t,x_init));

%% Solve trajectory optimization
% Tag
tag = [pwd, '/', name];
tolerance = [1e-6,1e-8];
source = pwd;
for n = 1:length(tolerance)
    % Change the folder
    target = ['Tolerance_',sprintf('%3.0E',tolerance(n))];
    mkdir(target);
    cd(target);
    % Run the optimization
    snoptOptions = setSolverOptions(snoptOptions,'MajorOptimalityTolerance',tolerance(n),'MajorFeasibilityTolerance',tolerance(n));
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

