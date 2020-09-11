function animateMultiSolution(plant, solns, labels, name)

% Get the total time of the solution
t = solns{1}.xtraj.getBreaks();
T = t(end);
% Resample the time axis
nFrames = round(30*T);
ts = linspace(0,T,nFrames);
% Resample the configurationn trajectories
nQ = plant.getNumPositions();
nsolns = length(solns);
q = cell(1,nsolns);
for n = 1:nsolns
   x = solns{n}.xtraj.eval(ts);
   q{n} = x(1:nQ,:);
end
% Animate the trajectory
animationUtilities.visualizeMultiTrajectory(plant, q, name, labels);


end