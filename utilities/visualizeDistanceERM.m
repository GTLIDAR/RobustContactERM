function  visualizeDistanceERM(plant, soln, sigma)
%VISUALIZEDISTANCEERM Summary of this function goes here
%   Detailed explanation goes here

% Get the time axis
t = soln.xtraj.getBreaks();
% Now get the configuration trajectory
x = soln.xtraj.eval(t);
nX = size(x,1);
q = x(1:nX/2,:);
% And get the normal forces
f = soln.ltraj.eval(t);
f = f(1,:);
% Calculate the normal distances
phi = zeros(size(f));
for n = 1:size(f,2)
   phi(:,n) = plant.contactConstraints(q(:,n)); 
end
% Expand the fixed sigma value
sigma = sigma * ones(size(phi));
% Visualize the results
ermCostVisualizer(plant, q, t, f, phi, sigma, 'Normal Force','Normal Distance');
end




