% Set up the falling rod example
model = FallingRod(1, 0.5, 0.05, 0.002, 0.6, 0.0025);
X0 = [0.0, 1.0, pi/6, 0.0, 0.0, 4.0]';

% Set the solver option
model.contactSolver = 'SLCP_Logistic';

sigmas = [0.01,0.1,1,10];

% Create three new figures
figs(1:3)   = figure('Name','Configuration');        % Configuration profiles    
figs(2)     = figure('Name','Velocity');             % Velocity profiles
figs(3)     = figure('Name','Contact Forces');       % Force Profiles

for n = 1:length(sigmas)
    % Set the value for Sigma (the variance)
    model.sigma = sigmas(n);
    % Run a simulation to get the trajectory
    [t, x, f, r] = model.simulate(X0, 1.0);
    % Plot the trajectories
    figs = plotTrajectories(figs,t,x,f,['\sigma = ',num2str(sigmas(n))]);
    % Plot the residual analysis
    plotResiduals(t,f,r,['\sigma = ',num2str(sigmas(n))]);
    % Plot the other diagnostics
    diagnostics(model,x,f,t,['sigma = ',num2str(sigmas(n))]);
    % Create an animation
%     frames = model.animate(x);
%     %Save the animation
%     w = VideoWriter(['fallingRod_SLCPLogistic_Sigma=',num2str(sigmas(n)),'.avi']); %#ok<TNMLP>
%     open(w);
%     for k = 1:length(frames)
%         writeVideo(w,frames(k))
%     end
%     close(w);
end

function plotResiduals(t,f,r,name)
figure();
yyaxis left;
f = sum(abs(f),1);
plot(t,f,'-');
ylabel('Force Magnitude');
ylim([-1,1]*max(f));
yyaxis right;
plot(t,r,'-');
ylim([-1,1]*max(abs(r)));
ylabel('Complementarity Residual');
xlabel('Time (s)');
title(name);
end

function figs = plotTrajectories(figs, t, x, f,name)
    
   % Plot the positions
   figure(figs(1));
   labels = {'COM_X','COM_Y','\theta'};
   for k = 1:3
      subplot(3,1,k)
      plot(t,x(k,:),'DisplayName',name);
      hold on;
      ylabel(labels{k});
   end
   xlabel('Time (s)');
   % Plot the velocities
   figure(figs(2));
   labels = {'Horizontal','Vertical','Angular'};
   for k = 1:3
    subplot(3,1,k);
    plot(t,x(3+k,:),'DisplayName',name);
    hold on;
    ylabel(labels{k});
   end
   xlabel('Time (s)');
   % Plot the forces
   % First, resolve the forces
    R = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 -1 0 0;
         0 0 0 0 1 -1];
    fr = R*f;
    % Now, plot the forces
    figure(figs(3));
    labels = {'Normal1','Normal2','Friction1','Friction2'};
    for k = 1:4
        subplot(4,1,k);
        plot(t,fr(k,:),'DisplayName',name);
        hold on;
        ylabel(labels{k});
    end
    xlabel('Time (s)')
end

function diagnostics(model, x, f, t, name)

% Calculate the distances from the ground
d = zeros(2, size(x, 2));
for k = 1:size(x,2)
    [n,a] = model.contactNormal(x(1:3,k));
    d(:,k) = n*x(1:3,k) - a; 
end
% Plot the Distances and Normal Forces in the same plot
figure('Name',name);
titles = {'Lower End','Higher End'};
for k = 1:2
    subplot(2,1,k);
    yyaxis left;
    plot(t,f(k,:),'-');
    ylabel('Normal Force');
    ylim([-0.1,1]*max(abs(f(k,:))));
    yyaxis right;
    plot(t,d(k,:),'-');
    ylim([-0.1,1]*max(abs(d(k,:))));
    ylabel('Distance');
    xlabel('Time (s)');
    title(titles{k});
end
end