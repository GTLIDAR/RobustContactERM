function soln = optimizePlant(plant, guess, optimOptions, solverOptions)
%% OptimizePlant: General purpose script for optimizing trajectories using Drake's trajectory optimization classes
%   optimizePlant optimizes contact-implicit mechanical system problems.
%
%
%   Arguments:
%       plant: a subclass of Drake's Manipulator, the dynamic system plant
%       used for optimization
%       
%       init_traj: a structure with the fields:
%           t:      initial time samples
%           traj:   initial guess at the trajectories
%
%       optimOptions: a structure with the fields:
%           durations: vector specifying the allowed durations for the
%               problem
%           npoints: the number of knot points for optimization
%           options: optimizer options, specific to the optimzer chosen
%           
%

%   February 27, 2020
%   Luke Drnach

%% Create the output structure
soln = struct('t',[],'xtraj',[],'utraj',[],'ltraj',[],'z',[],'F',[],'info',[],...
    'infeasible',[],'elapsed',[]);

%% Create the optimization problem
% Get the necessary parameters (duration and number of grid points) from the input structure
duration = optimOptions.duration;
npoints = optimOptions.nPoints;
options = optimOptions.options;
% Create the problem
fprintf('\nCreating Optimization Problem\n');
prob = RobustContactImplicitTrajectoryOptimizer(plant, npoints, duration, options);

%% Add in State constraints
if isfield(optimOptions, 'stateConstraints') && ~isempty(optimOptions.stateConstraints)
    fprintf('Adding state constraints\n');
    for n = 1:length(optimOptions.stateConstraints)
        x = optimOptions.stateConstraints(n).value;
        k = optimOptions.stateConstraints(n).knotPoint;
        idx = optimOptions.stateConstraints(n).stateIdx;
        % Create the constriant
        prob = prob.addStateConstraint(ConstantConstraint(x(idx)), k, idx);
    end
end

%% Add in Objectives
% Add in the running cost
if isfield(optimOptions, 'runningCost') && ~isempty(optimOptions.runningCost)
   fprintf('Adding running cost\n');
   if ~isstruct(optimOptions.runningCost)
       prob = prob.addRunningCost(optimOptions.runningCost);
   else
       for n = 1:length(optimOptions.runningCost)
          prob = prob.addRunningCost(optimOptions.runningCost(n).cost, optimOptions.runningCost(n).name); 
       end
   end
end
% Add in the final cost
if isfield(optimOptions, 'finalCost') && ~isempty(optimOptions.finalCost)
   fprintf('Adding terminal cost\n');
    prob = prob.addFinalCost(optimOptions.finalCost); 
end

%% Set the options for the solver
% Set the solver to SNOPT
prob = prob.setSolver('snopt');
% Set the solver options for the problem
solverFields = {'MajorFeasibilityTolerance','MajorOptimalityTolerance','MinorFeasibilityTolerance',...
    'ScaleOption','IterationsLimit','ElasticWeight','MajorIterationsLimit','SuperbasicsLimit','print'};
solve_fields = fieldnames(solverOptions);
fprintf('Setting solver options\n');
for n = 1:length(solve_fields)
    if any(strcmpi(solve_fields{n}, solverFields))
        prob = prob.setSolverOptions('snopt',solve_fields{n}, solverOptions.(solve_fields{n}));
    else
        fprintf('%s is not an allowable field for SNOPT and is ignored\n',solve_fields{n});
    end
end

% Enable or disable the display function
if isfield(optimOptions, 'display') && optimOptions.display
    fprintf('Objective and Constraint Display Enabled\n');
    prob = prob.enableCostDisplay();
end

%% Solve the optimization problem
% Get the initial time and trajectory from the input
t_init = guess.t;
traj_init = guess.traj;

fprintf('Running Trajectory Optimization\n');
pause(0.1);
tic;
[soln.xtraj, soln.utraj, soln.ltraj, soln.z, soln.F, soln.info, soln.infeasible] = prob.solveTraj(t_init, traj_init);
soln.elapsed = toc;
fprintf('Elapsed time is %f seconds\n',soln.elapsed);
pause(0.1);
% Check for infeasible constraints
if ~isempty(soln.infeasible)
    fprintf('Infeasible constraints:\n');
    cstr = unique(soln.infeasible);
    for n = 1:length(cstr)
       fprintf(' %s\n',cstr{n}); 
    end
end
% Get and store the time vector
dt = soln.z(prob.h_inds);
soln.t = cumsum([0;dt]);
% Check for relaxation variables
if optimOptions.options.nlcc_mode == 5
   soln.relax = soln.z(prob.relax_inds); 
   fprintf('Relaxed NC Constraint Maximum: %8.2e\n',max(soln.relax));
end
% Clear the persistent variables, if display is on
end