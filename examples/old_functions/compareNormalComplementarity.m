function figs = compareNormalComplementarity(plant, reference, ermTraj)
%COMPARENORMALCOMPLEMENTARITY Summary of this function goes here
%   Detailed explanation goes here


%   Luke Drnach
%   April 21, 2020

% Get the nominal data
t = reference.xtraj.getBreaks();
xN = reference.xtraj.eval(t);
lN = reference.ltraj.eval(t);

% Assume 2D, so we can construct the number of contact points
numL = size(lN, 1);
numContact = numL/4;
S = zeros(2*numContact,numL);        % Select the forces from the contact variables

S(1:numContact:end, 1:4:end) = eye(numContact);
S(2:numContact:end, 2:4:end) = eye(numContact);
S(2:numContact:end, 3:4:end) = -eye(numContact);
% nominal force trajectories
fN = S*lN;

% Prepare for convergence analysis
qN = xN(1:3,:);
phiN = zeros(1,size(qN,2));
for k = 1:size(qN,2)
    phiN(k) = plant.contactConstraints(qN(:,k));
end

phiErr = zeros(1,length(ermTraj));
forceErr = zeros(1,length(ermTraj));
sigmas  = zeros(1,length(ermTraj));
violation = zeros(1,length(ermTraj)); 
% Create the figure
figs(1) = figure();
for n = 1:length(ermTraj)
    t = ermTraj(n).xtraj.getBreaks();
    x = ermTraj(n).xtraj.eval(t);
    f = ermTraj(n).ltraj.eval(t);
    
    % Calculate the normal distances
    phi = zeros(1,size(x,2));
    for k = 1:length(phi)
       phi(k) = plant.contactConstraints(x(1:3,k)); 
    end
    subplot(2,1,1);
    plot(t,phi,'LineWidth',1.5,'DisplayName',sprintf('\\sigma = %0.2e',ermTraj(n).sigma));
    ylabel('Normal Distance');
    xlabel('Time (s)');
    hold on;
    
    subplot(2,1,2);
    plot(t(1:end-1),f(1,1:end-1),'LineWidth',1.5);
    ylabel('Normal Force');
    xlabel('Time (s)');
    hold on;
    
    % Calculate the average pointwise error
    phiErr(n) = mean((phiN - phi).^2);
    forceErr(n) = mean((fN(1,:) - f(1,:)).^2); 
    sigmas(n) = ermTraj(n).sigma;
    
    % Calculate the complementarity violation
    violation(n) = max(phi.*f(1,:));
    
end
% Add in the reference trajectory
subplot(2,1,1);

% Normal distance
plot(t,phiN,'k-','DisplayName','reference');
legend show;
legend boxoff;
% Normal force
subplot(2,1,2);
plot(t(1:end-1), fN(1,1:end-1),'k-');

%% ERM Convergence

% Check for info fields
if isfield(ermTraj,'info')
    status = [ermTraj(:).info];
else
    status = ones(1,ntraj);
end

% Convergence plots
figs(2) = figure();
subplot(2,1,1);
loglog(sigmas, phiErr, 'ko-','MarkerFaceColor','k','LineWidth',1.5);
ylabel('Normal Distance');
title('Deviation of ERM from Complementarity Solution');
if any(status > 1)
   hold on;
   loglog(sigmas(status > 1), phiErr(status > 1), 'ro');
end

subplot(2,1,2);
loglog(sigmas, forceErr, 'ko-','MarkerFaceColor','k','LineWidth',1.5);
ylabel('Normal Force');
xlabel(sprintf('Uncertainty (\\sigma)'));
if any(status > 1)
   hold on;
   loglog(sigmas(status > 1), forceErr(status > 1), 'ro');
end


% Complementarity violation
ref = max(phiN.*fN(1,:));
violation = [ref, violation];
sigmas = [0,sigmas];
figure();
plot(sigmas, violation, 'ko-','MarkerFaceColor','k','LineWidth',1.5);
ylabel('Max Compplementarity Violation');
xlabel('Uncertainty (\\sigma)');

end

