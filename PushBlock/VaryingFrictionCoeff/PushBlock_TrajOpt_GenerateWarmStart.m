function PushBlock_TrajOpt_GenerateWarmStart()
%PushBlock_TrajOpt Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   December 17, 2019

% Add in the necessary paths
here = pwd;
cd ..
add_drake;
cd(here);
addpath(genpath('PATH_LCP'));

name = 'PushBlock_TrajOpt_Nominal';
% Create the plant model
dt = 0.01;
plant = Block();
plant.timestep = dt;
plant.terrain.friction_coeff = 0.5;
% Specify the initial and final conditions
x0 = [0,plant.height/2,0,0]';     
xf = [5,plant.height/2,0,0]';    

% Specify the final time of the problem
Tf = 1;     % Five seconds for this problem
    
t_init = 0:plant.timestep:Tf;
N = length(t_init);
% Generate the initial guess trajectory
x_init = zeros(length(x0), length(t_init));
for n = 1:length(x0)
   x_init(n,:) = linspace(x0(n), xf(n), N);
end
traj_init.x = PPTrajectory(foh(t_init,x_init));
% Solve the optimization
[xtraj, utraj,z,F,info] = optimizePushBlock(plant, t_init,x0, xf, traj_init);
% Save the results
save(name,'plant','xtraj','utraj','z','F','info','x0','xf','t_init');
   
%% Visualize the results

% Get the results
[t,x] = getPointsFromTrajectory(xtraj);
[~,u] = getPointsFromTrajectory(utraj);

% Recalculate the contact forces
f0 = plant.contactForce(x(1:2,1),x(3:4,1),u(:,1));
f = zeros(length(f0),length(t));
f(:,1) = f0;
for n = 2:length(t)
   f(:,n) = plant.contactForce(x(1:2,n),x(3:4,n),u(:,n));
end

% Plot the results
figure();
subplot(3,1,1);
plot(t,x(1,:));
ylabel('Position');
subplot(3,1,2);
plot(t,x(3,:));
ylabel('Velocity');
subplot(3,1,3);
plot(t,u);
ylabel('Control');
xlabel('Time (s)');

% Contact impulses
figure();
subplot(2,1,1);
plot(t,f(1,:));
ylabel('Normal');
title('Contact Forces');
subplot(2,1,2);
plot(t,f(2,:) - f(3,:));
ylabel('Tangential');
xlabel('Time (s)');
end