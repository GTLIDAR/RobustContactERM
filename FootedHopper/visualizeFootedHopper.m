function  visualizeFootedHopper(plant, soln, savename)
%VISUALIZEFOOTEDHOPPER Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   May 11, 2020

% Force selector matrix
S = [1 0 0 0 0 0 0 0;
     0 1 -1 0 0 0 0 0;
     0 0 0 0 1 0 0 0;
     0 0 0 0 0 1 -1 0];
% Get the trajectories
t = soln.t;
x = soln.xtraj.eval(t);
u = soln.utraj.eval(t);
l = soln.ltraj.eval(t);
% Resolve the forces
f = S*l;

%% Plot the base trajectory
fig = figure('name','Base');
subplot(2,1,1);
plot(t,x(1,:),'LineWidth',1.5,'DisplayName','Horizontal');
hold on;
plot(t,x(2,:),'LineWidth',1.5,'DisplayName','Vertical');
ylabel('Base Position');
legend show;
legend boxoff;
subplot(2,1,2);
plot(t,x(6,:),t,x(7,:),'LineWidth',1.5);
ylabel('Base Velocity');
xlabel('Time (s)');
if nargin == 3
    saveas(fig,'BaseTrajectory.fig');
end

%% Plot the joint trajectories
fig = figure('name','Joints');
subplot(3,1,1);
hold on;
plot(t,x(3,:),'LineWidth',1.5,'DisplayName','Hip');
plot(t,x(4,:),'LineWidth',1.5,'DisplayName','Knee');
plot(t,x(5,:),'LineWidth',1.5,'DisplayName','Ankle');
legend show;
legend boxoff;
ylabel('Joint Angles');
subplot(3,1,2);
hold on;
plot(t,x(8,:),t,x(9,:),t,x(10,:),'LineWidth',1.5);
ylabel('Joint Velocity');
subplot(3,1,3);
hold on;
plot(t,u(1,:),t,u(2,:),t,u(3,:),'LineWidth',1.5);
ylabel('Joint Torques');
xlabel('Time (s)');
if nargin == 3
   saveas(fig,'JointTrajectory.fig'); 
end

%% Plot the force trajectories
fig = figure('name','Forces');
subplot(2,1,1);
hold on;
plot(t,f(1,:),'LineWidth',1.5,'DisplayName','Toe');
plot(t,f(3,:),'LineWidth',1.5,'DisplayName','Heel');
ylabel('Normal Force');
legend show;
legend boxoff;
subplot(2,1,2);
plot(t,f(2,:),t,f(4,:),'LineWidth',1.5);
ylabel('Friction Force');
xlabel('Time (s)');
if nargin ==3
   saveas(fig, 'ForceTrajectory.fig'); 
end

%% Animate the trajectory
q = x(1:5,:);
if nargin == 3
    utilities.trajectoryAnimator(plant, q, [], [savename,'.avi']);
else
    utilities.trajectoryAnimator(plant, q);
end

end

