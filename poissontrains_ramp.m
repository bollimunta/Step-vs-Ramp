%% generate poisson spike trains with ramping fr 
function [s]=poissontrains_ramp(ntrls,a,b,T)
% INPUT:
% ntrls - number of trials
% a - initial firing rate (spikes/sec)
% b - final firing rate (spikes/sec)
% T - time period [seconds]
% OUTPUT
% s - spike train (ntrls*time ; 1ms bins) matrix containing zeros(no spike) and
%     ones(one spike) or twos(two spikes)

dt=0.001;
s=zeros(ntrls,T/dt);
for i= 1:T/dt
    r=a + (b-a)*i*dt/T;
    temp=rand(ntrls,1);
    s(:,i)=temp<=r*dt;
end
end