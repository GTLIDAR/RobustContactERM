function sim_results = openLoopSimulationStudy()

% Luke Drnach
% February 7, 2020

% Load in the test case and the test controls
testCase = load('testCase.mat');
testCase = testCase.testCase;

controls = load('testControls.mat');
controls = controls.testControls;

% Get the simulation variables
target = testCase.xf;
plant = testCase.plant;
friction = testCase.friction;
x0 = testCase.x0;

% Run a series of simulations
% Initialize the final error structure
posdiff = zeros(length(controls), length(friction));
veldiff = zeros(length(controls), length(friction));

sim_results(length(controls)*length(friction)) = struct('state',[],'time',[],'control',[],'source',[],'uncertainty',[],...
    'x0',[],'target',[],'friction',[],'terminalError',[]);
labels = cell(1,length(controls));
i = 1;
for n = 1:length(controls)
    % Get the current control
    %u = controls(n).utraj.eval(controls(n).t_init);
    u = controls(n).utraj;
    Tf = controls(n).t_init(end);
    for k = 1:length(friction)
       % Change the plant friction
       plant.terrain.friction_coeff = friction(k);
       % Simulate the control
       [t, x] = plant.simulate(Tf, x0, u);
       % Get the achieved final state
       achieved = x(:,end);
       % Calculate the final state error
       sim_results(i).terminalError = achieved - target;
       posdiff(n,k) = sim_results(i).terminalError(1);
       veldiff(n,k) = sim_results(i).terminalError(3);
       
       % Cache the simulation results
       sim_results(i).state = x;
       sim_results(i).time = t;
       sim_results(i).control = u;
       sim_results(i).source = controls(n).case;
       sim_results(i).uncertainty = controls(n).sigma;
       sim_results(i).x0 = x0;
       sim_results(i).target = target;
       sim_results(i).friction = friction(k);
       i = i+1;
    end
    switch controls(n).case
        case 'nominal'
            labels{n} = 'mean';
        case 'ERM'
            labels{n} = sprintf('\\sigma = %0.4f',controls(n).sigma);
        otherwise
            labels{n} = controls(n).case;
    end
end
save('openLoopSimResults.mat','sim_results');


% Make a figure of the errors
figure();
idx = ones(size(posdiff));
idx = cumsum(idx);
subplot(2,1,1);
scatter(idx(:),posdiff(:),'filled');
ylabel('Position Difference');
set(gca, 'XTick',1:length(controls), 'XTickLabels',labels,'XTickLabelRotation',30);
xlim([0, length(controls)+1]);

subplot(2,1,2)
scatter(idx(:), veldiff(:), 'filled');
ylabel('Velocity Difference');
set(gca, 'XTick',1:length(controls), 'XTickLabels',labels,'XTickLabelRotation',30);
xlim([0, length(controls)+1]);
end

