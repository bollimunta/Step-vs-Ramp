% test the performance on artificial data

[s]=poissontrains_ramp(100,0,20,0.3); % generate test data that is ramping in firing rate
global_ratio_data(s) % run model selection
