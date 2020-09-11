function solutionEvaluation(soln)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Plot Complementarity Relationships
% Get the force variables
lambda = soln.ltraj.eval(soln.t);
nContact = size(soln.cstrVal.NormalDistanceCompl1, 1);
nFriction = size(soln.cstrVal.TangentVelocityCompl1,1)/nContact;
% Create Force Indices
tangent_inds = 1:size(lambda,1);
skip = 2*nFriction;
normal_inds = 1:skip:length(tangent_inds);
gamma_inds = nContact + nFriction : skip : length(tangent_inds);
tangent_inds([normal_inds, gamma_inds]) = [];
% First plot the complementarity conditions 
for n = 1:nContact
    % Normal Force Complementarity plot
    figure('name','Contact Constraints');
    subplot(nFriction + 2,1,1);
    title(sprintf('Contact Point %d',n));
    yyaxis left;
    plot(soln.t(1:end-1), lambda(normal_inds(n),1:end-1),'LineWidth',1.5);
    ylabel('Normal Force');
    yyaxis right;
    plot(soln.t(1:end-1), soln.contactVals.normalDistance(n,:),'LineWidth',1.5);
    ylabel('Normal Distance');
    xlabel('Time (s)');
    alignZero(gca);
    % Tangent Force Complementarity plots
    fric_inds = (n-1)*nFriction+1:n*nFriction;
    for k = 1:nFriction
        subplot(nFriction+2,1,k+1);
        yyaxis left;
        plot(soln.t(1:end-1), lambda(tangent_inds(fric_inds(k)), 1:end-1), 'LineWidth',1.5);
        ylabel('Friction Force');
        yyaxis right;
        plot(soln.t(1:end-1), soln.contactVals.slidingVelocity(fric_inds(k), :), 'LineWidth',1.5);
        ylabel('Tangent Velocity');
        xlabel('Time (s)');
        alignZero(gca);
    end
    % Sliding Velocity / Friction Cone Complementarity Plot
    subplot(nFriction + 2,1,nFriction + 2);
    yyaxis left;
    plot(soln.t(1:end-1), lambda(gamma_inds(n), 1:end-1), 'LineWidth',1.5);
    ylabel('Sliding Velocity Slack');
    yyaxis right
    plot(soln.t(1:end-1), soln.contactVals.frictionCone(n,:),'LineWidth',1.5);
    ylabel('Friction Cone Defect');
    xlabel('Time (s)');
    alignZero(gca);
end
% Plot the constraints as seen by the optimizer
figure('name','Enforced Constraints');
subplot(nFriction + 2, 1, 1);
plot(soln.t(1:end-1), soln.cstrVal.NormalDistanceNonNeg1,'LineWidth',1.5);
ylabel('Normal Distance > 0');
for k = 1:nFriction
    subplot(nFriction + 2, 1, 1 + k)
    plot(soln.t(1:end-1), soln.cstrVal.TangentVelocityNonneg1(k:nFriction:end,:),'LineWidth',1.5);
    ylabel(sprintf('Tangent Velocity %d > 0', k));
end
subplot(nFriction+2,1,nFriction+2)
plot(soln.t(1:end-1), soln.cstrVal.FrictionConeNonneg1,'LineWidth',1.5);
ylabel('Friction Cone > 0');
xlabel('Time (s)');

% Plot the defects in another figure
figure('name','Complementarity Defects');
subplot(3,1,1);
plot(soln.t(1:end-1), soln.cstrVal.NormalDistanceCompl1,'LineWidth',1.5);
ylabel('Normal Distance Residual');
subplot(3,1,2);
plot(soln.t(1:end-1), soln.cstrVal.TangentVelocityCompl1,'LineWidth',1.5);
ylabel('Tangent Velocity Residual');
subplot(3,1,3)
plot(soln.t(1:end-1), soln.cstrVal.FrictionConeCompl1,'LineWidth',1.5);
ylabel('Friction Cone Residual');
xlabel('Time (s)');

%% Plot Dynamic Constraints
figure('Name','Dynamics Constraint Violations');
nQ = size(soln.cstrVal.DynamicConstraints1, 1)/2;
subplot(2,1,1);
plot(soln.t(1:end-1), soln.cstrVal.DynamicConstraints1(1:nQ,:), 'LineWidth',1.5);
ylabel({'Dynamic Defects','Configuration'});
subplot(2,1,2);
plot(soln.t(1:end-1), soln.cstrVal.DynamicConstraints1(nQ+1:end,:),'LineWidth',1.5);
ylabel({'Dynamic Defects','Velocity'});
xlabel('Time (s)');

%% Plot Joint Limit Forces
if isfield(soln, 'jltraj')
    figure('Name','Joint Limit Force Evaluation');
    jl = soln.jltraj.eval(soln.t);
    nJL = size(jl, 1)/2;
    for n = 1:nJL
       subplot(nJL, 1, n);
       plot(soln.t, jl(n,:), soln.t, jl(nJL+n,:),'LineWidth',1.5);
       ylabel({'Joint Limit',sprintf('Force %d',n)}); 
    end
end

%% Plot Cost Over Time
if isfield(soln,'costs')
    figure('name','Running Costs');
    costNames = fieldnames(soln.costs);
    nCosts = length(costNames);
    for n = 1:nCosts
        subplot(nCosts, 1, n);
        plot(soln.t(1:end-1), soln.costs.(costNames{n}),'LineWidth',1.5);
        ylabel(costNames{n});
    end
    xlabel('Time (s)');
end
end


function ax = alignZero(ax)
yyaxis(ax, 'left');
left_lims = ylim;
yyaxis(ax,'right');
right_lims = ylim;

l_range = range(left_lims);
r_range = range(right_lims);

l_percent = left_lims/l_range;
r_percent = right_lims/r_range;

new_percent = [min(l_percent(1), r_percent(1)), max(l_percent(2), r_percent(2))];

yyaxis(ax, 'left')
ylim(new_percent * l_range);
yyaxis(ax, 'right');
ylim(new_percent * r_range);



end

