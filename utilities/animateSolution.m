function animateSolution(plant, soln, name)
%ANIMATESOLUTION 

%   Resamples solution trajectory at 30Hz and makes the corresponding
%   trajectory animation
%

% Luke Drnach
% August 31, 2020

% Get the total time of the solution
t = soln.xtraj.getBreaks();
T = t(end);

% resample at 30Hz
Nframes = round(T*30);
ts = linspace(0,T,Nframes);

% Get the configurations
x = soln.xtraj.eval(ts);
nQ = plant.getNumPositions();
q = x(1:nQ,:);

% Make the animation
if nargin == 3
    animationUtilities.visualizeTrajectory(plant, q, name);
else
    animationUtilities.visualizeTrajectory(plant, q);
end

end

