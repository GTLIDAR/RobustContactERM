function compareControlEnergy(nominal, ermTraj)
% compareControlEnergy: Compares the control energy used by the ERM
% Trajectories to the control energy used by the nominal model

%   Luke Drnach
%   February 19, 2020

nTraj = length(ermTraj);
% Get the nominal data
t = nominal.xtraj.getBreaks();
dt = diff(t);

uN = nominal.utraj.eval(t);

% Initialize arrays
ctrlNRG = zeros(1,nTraj+1);
ctrlNRG(1) = energy(uN(:,1:end-1), dt);

label = cell(1,nTraj+1);
label{1} = 'nominal';

% Loop over all the data
for n = 1:nTraj
    % Get the trajectories from the current run
    u = ermTraj(n).utraj.eval(t);     % Controls
    % calculate the errors
    ctrlNRG(n+1) = energy(u(:,1:end-1), dt);
    label{n+1} = sprintf('\\sigma = %0.4f',ermTraj(n).sigma);
end

% Make the final plot
figure();
bar(ctrlNRG);
ylabel('Control Energy');
set(gca, 'XTick',1:nTraj+1, 'XTickLabel',label, 'XTickLabelRotation',30);
end

function nrg = energy(u, dt)
U = diag(diag(u'*u));
nrg = sum(U*dt(:));
end