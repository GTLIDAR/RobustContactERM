function  FootedHopper_Flat_Feasible_Slack_VariableTimesteps()
%FOOTEDHOPPER_JOINTLIMIT_FEASIBLE Summary of this function goes here
%   Detailed explanation goes here
%% Set-up for the test
% Name the test
stk = dbstack;
name = stk.name;
tag = [pwd,'/',name];
% Get the low-level options for SNOPT
snoptOptions = setSolverOptions([], 'IterationsLimit',50000','MajorIterationsLimit',5000,'SuperbasicsLimit',1500,'ScaleOption',1);
% Get the high-level options for the optimization
optimOptions = setOptimizerOptions([], 'nPoints',101, 'integration_method',RobustContactImplicitTrajectoryOptimization.BACKWARD_EULER, ...
    'nlcc_mode',2,'duration',[1,5],'distanceScaling',1, 'display',true,'time_option',2);

%% Create the plant model
plant = FootedHopper();
% Add Joint Limits;
jl_min = [-inf, -inf, -90, -120, 30]' * pi / 180;
jl_max = [inf, inf, 90, 60, 150]' * pi / 180;
plant = plant.setJointLimits(jl_min, jl_max);
plant.lengths(3) = 1/3;
plant.masses(3) = 1/3;
plant.inertias(3) = plant.masses(3)*plant.lengths(3)^2 /12;

%% Change the terrain to a step terrain
plant.terrain = FlatTerrain();

%% Set the boundary conditions
% Enpoint boundary conditions
x0 = [0, 1.5, (1-plant.heelFraction)*plant.lengths(3),0]';
xf = [4, 1.5, 4 + (1-plant.heelFraction)*plant.lengths(3), 0]';
% Calculate the initial and final states
q0 = plant.inverseKinematics(x0(1:2), x0(3:4),0);
qf = plant.inverseKinematics(xf(1:2), xf(3:4),0);
% Set the desired velocities to zero
x0 = [q0(:); zeros(5, 1)];
xf = [qf(:); zeros(5, 1)];
% Set the indices for partial state constraints
x0_idx = 1:10;
xf_idx = 1:10;
% Create constraints
optimOptions.stateConstraints(1) = setStateConstraint(x0, 1, x0_idx);
optimOptions.stateConstraints(2) = setStateConstraint(xf, optimOptions.nPoints, xf_idx);

%% Generate an initial guess for the trajectory
% Use linear interpolation between the initial and final states
t_init = linspace(0, optimOptions.duration(2), optimOptions.nPoints);
% linear interpolated warmstart
x_init = zeros(10, optimOptions.nPoints);
% for n = 1:10
%    x_init(n,:) = linspace(x0(n), xf(n), optimOptions.nPoints); 
% end
% Store in a structure
guess.t = t_init;
guess.traj.x = PPTrajectory(foh(t_init, x_init));

slacks = [10,1,0.1,0.01,0.001,0.0001,0];
errorIdx = false(size(slacks));
parent = pwd;
%% Solve the problem
for n = 1:length(slacks)
    % Change the folder
    target = ['Slack_',sprintf('%3.0E',slacks(n))];
    mkdir(target);
    cd(target);
    % Solve the problem
    optimOptions = setOptimizerOptions(optimOptions, 'compl_slack',slacks(n));
    soln  = optimizePlant(plant, guess, optimOptions, snoptOptions);
    
    %% Save the results and visualize
    save([name, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');
    
    % Print a report
    printReport(plant, soln, snoptOptions, optimOptions,'zeros');
    
    % Visualizer
    try
        visualizeFootedHopper(plant, soln, name);
    catch
        errorIdx(n) = true;
    end
    % Save the iteration figures
    savefig(figure(1),'CostandConstraintIterations.fig');
    % Clear all figures
    close all;
    % Change the initial guess
    guess.traj.x = soln.xtraj;
    guess.traj.u = soln.utraj;
    guess.traj.lambda = soln.ltraj;
    guess.traj.slacks = soln.slacks;
    guess.traj.jltraj = soln.jltraj;
    
    % Move to the parent directory
    cd(parent);
    
end

if any(errorIdx)
   for n = 1:length(errorIdx)
       if errorIdx(n)
           fprintf('Animation creation errors for slack %3.0E\n',slacks(n));
       end
   end
end
end