function fallingBlockSimulation()
%FALLINGBLOCKSIMULATION: Runs a simulation of a falling block

%   Luke Drnach
%   December 17, 2019

% Create the block model
plant = Block();
dt = plant.timestep;
% Initial condition
x0 = [0, 5, 1, 0]';  % Block starts 5 units off the terrain, with only horizontal velocity

% Final time 
Tf = 5;     % Simulate for 5 seconds
plant = plant.setupLCPCache(0:dt:Tf);

% Run the simulation
[t,x] = plant.simulate(Tf, x0);

% Create an animation of the results
q = x(1:2,:);
utilities.trajectoryAnimator(plant, q, [], 'FallingBlockSim.avi');

f = plant.lcpCache.data.force;

% Plot the trajectory
figure();
subplot(2,1,1);
plot(t,x(1,:));
title('Configuration');
ylabel('Horizontal Position');

subplot(2,1,2);
plot(t,x(2,:));
ylabel('Vertical Position');
xlabel('Time (s)');

% VELOCITY
figure();
subplot(2,1,1);
plot(t,x(3,:));
title('Velocity');
ylabel('Horizontal Velocity');

subplot(2,1,2);
plot(t,x(4,:));
ylabel('Vertical Velocity');
xlabel('Time (s)');

% FORCE
figure();
subplot(2,1,1);
plot(t,f(1,:));
title('Contact Force');
ylabel('Normal Force');

subplot(2,1,2);
plot(t,f(2,:) - f(3,:));
ylabel('Tangential Force');
xlabel('Time (s)');
end

