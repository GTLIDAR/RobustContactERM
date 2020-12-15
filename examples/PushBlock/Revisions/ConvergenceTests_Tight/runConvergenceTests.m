source = pwd;
cd Nominal
PushBlock_Nominal();
cd(source);
cd ERM_Warmstart
PushBlock_ERM_Warmstart();
cd(source);
cd ERM_NoWarmstart
PushBlock_ERM_NoWarmstart();
cd(source);
cd ERM_WeightTest
PushBlock_ERM_WeightTest();
cd(source);