function  DistanceScaled_Joint1Penalty_RobustComparison()
%DISTANCESCALED_ROBUSTCOMPARISON Summary of this function goes here
%   Detailed explanation goes here

% Load the ERM Solutions and the 'Nearest Optimal' Solution;
ermSolns = load([pwd,'/ContactCart_DistanceScaled_Joint1Penalty_ERM.mat']);
model = ermSolns.plant;

%% Make the convergence plots

% First compare the ERM Solns to the WarmStart
warmstart = convertWarmstart(ermSolns.guess);
ermTraj = ermSolns.soln;

compareRobustTrajectories(warmstart, ermTraj);

%% Compare the normal distance and normal forces
compareNormalComplementarity(model, warmstart, ermTraj);

%% Make an animation of all the ERM solutions with the non-ERM solutions together

target = ermSolns.optimOptions.stateConstraints(2).value(1:3);

q = cell(1, length(ermTraj)+1);
labels = cell(1,length(ermTraj) + 1);
for n = 1:length(ermTraj)
   x = ermTraj(n).xtraj.eval(ermTraj(n).t);
   q{n} = x(1:3,:);
   labels{n} = sprintf('\\sigma = %0.2E',ermTraj(n).sigma); 
end

% For the warm start
x = warmstart.xtraj.eval(warmstart.xtraj.getBreaks());
q{end} = x(1:3,:);
labels{end} = 'reference';
savename = 'DistanceScaledERM_Joint1Penalty.avi';
animationUtilities.multiTrajectoryAnimator(model, q, target, savename, labels);
end

function warmstart = convertWarmstart(guess)
    warmstart = struct('xtraj',guess.traj.x,'utraj',guess.traj.u,'ltraj',guess.traj.lambda);
end
