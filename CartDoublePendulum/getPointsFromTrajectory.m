function [t, x] = getPointsFromTrajectory(traj)
%GETPOINTSFROMTRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

t = traj.getBreaks();

x0 = traj.eval(t(1));

x = zeros(length(x0), length(t));
x(:,1) = x0;

for n = 2:length(t)
    x(:,n) = traj.eval(t(n));
end
end

