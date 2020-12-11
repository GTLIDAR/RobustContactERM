function ContactCart_DistanceScaled_Feasible()
%UNTITLED20 Summary of this function goes here
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
    'nlcc_mode',2,'duration',[1,1],'distanceScaling',10);

%% Create the plant model
plant = ContactDrivenCart();
plant.cartHeight = 1.5;         % Set the cart up higher
plant.inertias = [1/12, 1/12];

%% Set the boundary conditions
% Enpoint boundary conditions
x0 = [0,0]';
xf = [5, 0]';
% Calculate the initial and final states
q0 = plant.inverseKinematics(x0);
qf = plant.inverseKinematics(xf);
% Set the desired velocities to zero
x0 = [q0; zeros(3, 1)];
xf = [qf; zeros(3, 1)];
% Set the indices for partial state constraints
x0_idx = 1:6;
xf_idx = 1:6;
% Create constraints
optimOptions.stateConstraints(1) = setStateConstraint(x0, 1, x0_idx);
optimOptions.stateConstraints(2) = setStateConstraint(xf, optimOptions.nPoints, xf_idx);

%% Generate an initial guess for the trajectory
% Use linear interpolation between the initial and final states
t_init = linspace(0, optimOptions.duration(2), optimOptions.nPoints);
x_init = contactCartInitialGuess(x0, xf, t_init, 'linear');
% Store in a structure
guess.t = t_init;
guess.traj.x = PPTrajectory(foh(t_init, x_init));

%% Solve the trajectory optimization problem
soln = optimizePlant(plant, guess, optimOptions, snoptOptions);

%% Save and visualize the results.
save([name, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');

% Print a report
printReport(plant, soln, snoptOptions, optimOptions, 'linearInterpolation');

% Visualizer
visualizeCartTrajectories(plant, soln, name);

end

