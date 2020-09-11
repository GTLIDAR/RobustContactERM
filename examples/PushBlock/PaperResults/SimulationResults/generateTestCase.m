function testCase = generateTestCase()
%% generateTestCase: Generates the Test Conditions for running the Simulation Tests
%
%   Luke Drnach
%   February 10, 2020

% Load the nominal trajectory
nominal = load('../WarmStarts/PushBlock_Robust_NCP_NominalRefined.mat');

% Create a test case from the nominal model, using the plant, initial
% conditions, and final conditions
testCase.plant = nominal.plant;
testCase.x0 = nominal.x0;
testCase.xf = nominal.xf;

% Now create multiple values of friction to test
mu_min = 0.3;
mu_max = 0.7;
N_mu = 10;

testCase.friction = linspace(mu_min, mu_max, N_mu);

if nargout == 0
   save('testCase.mat','testCase'); 
end


end