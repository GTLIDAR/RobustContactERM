function distanceERMFigureMaker()
%DISTANCEERMFIGUREMAKER Summary of this function goes here
%   Detailed explanation goes here
S = [1, 0, 0, 0;
     0, 1,-1, 0];
% Load the ERM Solutions and the 'Nearest Optimal' Solution;
ermSolns = load([pwd,'/ContactCart_DistanceScaled_Joint1Penalty_ERM.mat']);
model = ermSolns.plant;
model.cart_height =0.5;
model.cart_width = 0.5;
% First Compare the ERM Trajectories to the WarmStart
warmstart = convertWarmstart(ermSolns.guess);
ermTraj = ermSolns.soln;

numErm = length(ermTraj);

% Get the Warmstart Trajectories
t = warmstart.xtraj.getBreaks();
xN = warmstart.xtraj.eval(t);
fN = S*warmstart.ltraj.eval(t);
uN = warmstart.utraj.eval(t);

phiN = getDistance(model, xN);

% Initialize error arrays
stateErr = zeros(1,numErm);
controlErr = zeros(1,numErm);
forceErr = zeros(1,numErm);
distanceErr = zeros(1,numErm);
normalErr = zeros(1,numErm);
sigmas = zeros(1,numErm);

avgDistance = zeros(1,numErm);
maxDistance = zeros(1,numErm);
feasibility = zeros(1,numErm);

% Loop over the ERM Solutions and plot them (the normal distance / normal
% force)
q = cell(1,numErm+1);
label = cell(1,numErm+1);

fig = figure();
for k = 1:numErm
    % Get the ERM Solution
    x = ermTraj(k).xtraj.eval(t);
    u = ermTraj(k).utraj.eval(t);
    f = S*ermTraj(k).ltraj.eval(t);
    % Get the normal distance trajectory 
    phi = getDistance(model, x);
    sigmas(k) = ermTraj(k).sigma;
    % Plot normal distance and normal force
    label{k} = sprintf('\\sigma = %0.3E',sigmas(k));
    figure(fig);
    subplot(2,1,1);
    plot(t,phi,'LineWidth',1.5,'DisplayName',label{k});
    hold on;
    subplot(2,1,2);
    plot(t,f(1,:),'LineWidth',1.5);
    hold on;
    % Calculate errors
    stateErr(k) = mean(sum((x - xN).^2, 1));
    controlErr(k) = mean(sum((u-uN).^2, 1));
    forceErr(k) = mean(sum((f - fN).^2, 1));
    distanceErr(k) = mean(sum((phi-phiN).^2,1));
    normalErr(k) = mean(sum((f(1,:) -fN(1,:)).^2, 1));
    % Save the configuration trajectory for frame animation
    q{k} = x(1:3,:);
    % Calculate some values for a table
    avgDistance(k) = mean(phi);
    maxDistance(k) = max(phi);
    feasibility(k) = max(phi.*f(1,:));
end
q{end} = xN(1:3,:);

% Add the WarmStart Trajectory
figure(fig);
subplot(2,1,1);
plot(t,phiN,'LineWidth',1.5,'DisplayName','reference');
ylabel('Normal Distance');
legend show;
legend boxoff;
subplot(2,1,2);
plot(t,fN(1,:),'LineWidth',1.5);
ylabel('Normal Force');
xlabel('Time');

% Create a figure for the convergence plots
figure();
subplot(2,1,1);
ln = loglog(sigmas, stateErr, 'o-','LineWidth',1.5,'DisplayName','State');
ln.MarkerFaceColor = ln.Color;
hold on;
ln = loglog(sigmas, controlErr, 'o-','LineWidth',1.5,'DisplayName','Control');
ln.MarkerFaceColor = ln.Color;
ln = loglog(sigmas, forceErr,'o-','LineWidth',1.5,'DisplayName','Force');
ln.MarkerFaceColor = ln.Color;
xlabel('Uncertainty (\sigma)');
ylabel('Mean-Square Difference');
legend show;
legend boxoff;

subplot(2,1,2);
ln = loglog(sigmas, distanceErr, 'o-','LineWidth',1.5,'DisplayName','Distance');
ln.MarkerFaceColor = ln.Color;
hold on;
ln = loglog(sigmas, normalErr,'o-','LineWidth',1.5,'DisplayName','Normal Force');
ln.MarkerFaceColor = ln.Color;
xlabel('Uncertainty (\sigma)');
ylabel('Mean-Square Difference');
legend show;
legend boxoff;

% Make and print a table of distance values
maxDistance = [max(phiN), maxDistance];
avgDistance = [mean(phiN), avgDistance];
feasibility = [max(phiN.*fN(1,:)), feasibility];
sigmas = [0, sigmas];

T = table(sigmas', maxDistance', avgDistance', feasibility','VariableNames',{'Sigma','MaxDistance','MeanDistance','Feasibility'});
disp(T);
writetable(T,'ErmDistanceData.csv');

% Produce a figure of some keyframes from the trajectory
label{end} = 'Reference';
% target = ermSolns.optimOptions.stateConstraints(2).value(1:3);
animationUtilities.visualizeFrames(model, q(end-3:end), 5, label(end-3:end));


end

function warmstart = convertWarmstart(guess)
    warmstart = struct('xtraj',guess.traj.x,'utraj',guess.traj.u,'ltraj',guess.traj.lambda);
end

function phi = getDistance(plant, x)
phi = zeros(1,size(x,2));
for n = 1:size(x,2)
   phi(n) = plant.contactConstraints(x(1:3,n)); 
end
end
