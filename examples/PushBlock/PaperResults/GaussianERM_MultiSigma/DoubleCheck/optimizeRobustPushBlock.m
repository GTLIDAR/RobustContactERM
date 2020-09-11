function [xtraj, utraj,ltraj,z,F,info,costVals] = optimizeRobustPushBlock(plant, t_init,x0, xf, traj_init, options)

% Get the final time and number of knot points
Tf = t_init(end);
N = length(t_init);
% Create a trajectory optimization problem
% Create the problem to solve
prob = RobustContactImplicitTrajectoryOptimizer(plant, N, Tf, options);
% Add the running cost
prob = prob.addRunningCost(@(t, x, u)cost2(t, x, u, xf));
% Add the initial and final value constraints
prob = prob.addStateConstraint(ConstantConstraint(x0), 1);
prob = prob.addStateConstraint(ConstantConstraint(xf), N);
% Set options for the solver
prob = prob.setSolver('snopt');
prob = prob.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-8);     %Default: 1e-6
prob = prob.setSolverOptions('snopt','MajorOptimalityTolerance',1e-8);      %Default: 1e-6
prob = prob.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-6);     %Default: 1e-6
prob = prob.setSolverOptions('snopt','ScaleOption',1);                      %Default: 1
prob = prob.setSolverOptions('snopt','IterationsLimit',50000);              %Default: 10,000
prob = prob.setSolverOptions('snopt','ElasticWeight',10^4);                 %Default: 10^4
prob = prob.setSolverOptions('snopt','MajorIterationsLimit',5000);          %Default: 1000
prob = prob.setSolverOptions('snopt','SuperbasicsLimit',1500);

if nargout == 7
    %prob = prob.enableDisplayFunction();
    initVals = struct('running',[],'dynamic',[],'distance',[],'friction',[],'sliding',[]);
    % Get the initial costs
    [initVals.running, initVals.dynamic, initVals.distance, initVals.friction, initVals.sliding] = ...
        prob.calculateCosts(diff(t_init), traj_init.x.eval(t_init), traj_init.u.eval(t_init), traj_init.lambda.eval(t_init));
end

% Solve the problem
disp('Running Trajectory Optimization');
pause(0.2);
tic;
[xtraj, utraj, ltraj, z, F, info] = prob.solveTraj(t_init, traj_init);
toc
pause(0.2);

if nargout == 7
    % Get the final costs
    finalVals = struct('running',[],'dynamic',[],'distance',[],'friction',[],'sliding',[]);
    [finalVals.running, finalVals.dynamic, finalVals.distance, finalVals.friction, finalVals.sliding] = ...
        prob.calculateCosts(diff(t_init), xtraj.eval(t_init), utraj.eval(t_init), ltraj.eval(t_init));
    costVals = struct('initial',initVals, 'final',finalVals);
end
end
% 
% function [g, dg] = cost(dt, x, u)
% % Running cost
% R = eye(numel(u)); % Control weights
% Q = eye(numel(x)); % State weights
% g = 1/2 * (x'*Q*x + u'*R*u)';
% % The differential cost
% dg = [0, x'*Q, u'*R];
% 
% end

function [g, dg] = cost2(~, x, u, xf)
% Running cost
R = 100*eye(numel(u)); % Control weights
Q = eye(numel(x)); % State weights

g = 1/2 * ( (x - xf)'*Q*(x - xf) + u' * R * u);
% The differential cost
dg = [zeros(1),(x - xf)'*Q,u'*R];
end

% function [h, dh] = terminalCost(t, x)
% % Terminal Cost
% h = t;
% % Differential Terminal Cost
% dh = [1, zeros(1,size(x,1))];
% end