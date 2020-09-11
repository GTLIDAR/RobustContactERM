function  continuehopperopt(filename )
old = load(filename);
% Pull the pland and options
plant = old.plant;
optimOptions = old.optimOptions;
snoptOptions = old.snoptOptions;
tag = old.tag;

% Pull the old solution to use as a guess
guess.t = old.soln.xtraj.getBreaks();
guess.traj.x = old.soln.xtraj;
guess.traj.u = old.soln.utraj;
if isfield(old.soln, 'ltraj')
    guess.traj.lambda = old.soln.ltraj;
end
if isfield(old.soln, 'jltraj')
    guess.traj.jltraj = old.soln.jltraj;
end
if isfield(old.soln,'slacks')
    guess.traj.slacks = old.soln.slacks;
end

% reset the running costs
xf = optimOptions.stateConstraints(2).value;
%optimOptions.runningCost = @(t, x, u)running_cost(t, x, u, xf, plant);
optimOptions.runningCost(1).cost = @(t,x,u)controlCost(t,x,u);
optimOptions.runningCost(1).name = 'ControlCost';
optimOptions.runningCost(2).cost = @(t,x,u)stateCost(t,x,u,xf);
optimOptions.runningCost(2).name = 'StateCost';
% Continue the optimization with the old run settings
soln  = optimizePlant(plant, guess, optimOptions, snoptOptions);
mkdir('continued');
cd continued
% Print the report
printReport(plant, soln, snoptOptions, optimOptions,filename,'continued_report.txt');
% Save the new results
parts = split(filename,'.');
new = [parts{1},'_continued'];

% Pull and save the results
save([new, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');
% Visualize the results
visualizeFootedHopper(plant, soln, new);

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
