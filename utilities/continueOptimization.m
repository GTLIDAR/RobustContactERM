function  continueOptimization(filename)
%CONTINUEOPTIMIZATION Summary of this function goes here
%   Detailed explanation goes here

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


% Continue the optimization with the old run settings
soln  = optimizePlant(plant, guess, optimOptions, snoptOptions);

% Print the report
printReport(plant, soln, snoptOptions, optimOptions,filename,'continued_report.txt');
% Save the new results
parts = split(filename,'.');
new = [parts{1},'_continued'];

% Pull and save the results
save([new, '.mat'], 'plant', 'guess', 'optimOptions', 'snoptOptions', 'soln', 'tag');



end
