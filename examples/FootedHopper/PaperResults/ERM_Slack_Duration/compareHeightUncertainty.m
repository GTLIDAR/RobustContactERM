% Load the data
labels = {'Sigma_5E-02','Sigma_9E-02','Sigma_3E-01','Sigma_5E-01'};
erm = cell(1,length(labels));
for n = 1:4
    data = load(['3Second_Combined_HeightSweep/ERM_HeightUncertainty_',labels{n},'/FootedHopper_Flat_3Second_Combined_HeightSweep.mat']);
    erm{n} = data.soln;
end
labels = {'\sigma = 0.05', '\sigma = 0.09', '\sigma = 0.30','\sigma = 0.50'};
% Load the reference solution
ref.xtraj = data.guess.traj.x;
% Compare the solutions
compareFootedHopperERM(data.plant, ref, erm, labels, '3Second_Combined_HeightSweep_Comparison');
