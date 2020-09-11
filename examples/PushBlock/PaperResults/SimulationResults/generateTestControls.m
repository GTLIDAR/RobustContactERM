function testControls = generateTestControls()
% GENERATETESTCONTROLS: Script to collect data for use in simulation tests.

%   Luke Drnach
%   February 10, 2020


% First load the nominal data
nominal = load('../WarmStarts/PushBlock_Robust_NCP_NominalRefined.mat');
testControls = addControl([], nominal,'nominal');

% Then load the 'worst case' baseline data
baseline = load('../Baseline_MultiFriction/PushBlock_MultiFrictionTraj_Robust_Baseline.mat');
baseline = baseline.trajData;
mu = [baseline(:).friction];
% Add the 'min' friction 
testControls = addControl(testControls, baseline(mu == 0.3), 'min');
% Add the 'max' friction
testControls = addControl(testControls, baseline(mu == 0.7), 'max');

% Finally, load the ERM data
uncertainty = load('../GaussianERM_MultiSigma/PushBlock_MultiSigma_ERM_Gaussian.mat');
uncertainty = uncertainty.trajData;
% sigma = [uncertainty(:).sigma];
% testVals = [0.001, 0.01, 0.1, 1.0];
for n = 1:length(uncertainty)
   testControls = addControl(testControls, uncertainty(n), 'ERM');
end

if nargout == 0
   save('testControls.mat','testControls'); 
end

end

function testControls = addControl(testControls, example, label)
%% addControl: Helper function for adding the relevant information to the testControls structure


sample.t_init = example.utraj.getBreaks();
% Add the controls
sample.utraj = example.utraj.eval(sample.t_init);
% Add the state trajectory
sample.xtraj = example.xtraj.eval(sample.t_init);
% Add the contact forces
sample.ftraj = example.ltraj.eval(sample.t_init);
% Add the label
sample.case = label;
% Add the friction coefficient
sample.friction = example.plant.terrain.friction_coeff;
% If there is one, add the uncertainty
if isfield(example, 'sigma')
    sample.sigma = example.sigma;
else
    sample.sigma = 0;
end

% Return the test controls
if isempty(testControls)
    testControls  =sample;
else
    testControls = [testControls; sample];
end
end
