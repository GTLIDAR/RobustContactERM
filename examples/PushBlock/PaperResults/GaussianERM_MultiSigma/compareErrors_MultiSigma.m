function compareErrors_MultiSigma()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%   Luke Drnach
%   February 7, 2020

% Force selector matrix
S = [1,0,0,0;
     0,1,-1,0];

% Load the nominal trajectory
nominal = load('../WarmStarts/PushBlock_Robust_NCP_NominalRefined.mat');

% Get the nominal trajectory
t = nominal.t_init;
x_nominal = nominal.xtraj.eval(t);
u_nominal = nominal.utraj.eval(t);
f_nominal = S*nominal.ltraj.eval(t);

% Load the ERM trajectories
ermTraj = load('PushBlock_MultiSigma_ERM_Gaussian.mat');
ermTraj = ermTraj.trajData;
% Store how much the trajectories diverge from the nominal
x_errors = zeros(1,length(ermTraj));
u_errors = zeros(1,length(ermTraj));
f_errors = zeros(1,length(ermTraj));
sigmas = zeros(1,length(ermTraj));

for n = 1:length(ermTraj)
   % Get the individual ERM trajectories
   x = ermTraj(n).xtraj.eval(t);
   u = ermTraj(n).utraj.eval(t);
   f = S * ermTraj(n).ltraj.eval(t);
   % Now calculate the divergences
   % State divergence
   xdiv = sum((x - x_nominal).^2, 1);
   x_errors(n) = sum(xdiv)./numel(xdiv);
   % Control divergence
   udiv = sum((u - u_nominal).^2, 1);
   u_errors(n) = sum(udiv)./numel(udiv);
   % Force divergence
   fdiv = sum((f - f_nominal).^2, 1);
   f_errors(n) = sum(fdiv)./numel(fdiv);
   % Standard Deviation
   sigmas(n) = ermTraj(n).sigma;
end

% Plot all three divergences
figure('name','Trajectory Divergence');
subplot(3,1,1);
loglog(sigmas, x_errors, 'ko-','MarkerFaceColor','k');
ylabel('State');
title('Trajectory Difference from Nominal');
subplot(3,1,2);
loglog(sigmas, u_errors, 'ko-','MarkerFaceColor','k');
ylabel('Control');
subplot(3,1,3);
loglog(sigmas, f_errors, 'ko-','MarkerFaceColor','k');
ylabel('Force');
xlabel('Friction Uncertainty');

end

