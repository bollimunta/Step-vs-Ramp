%% generate poisson spike trains with stepping fr 
function [s,n]=poissontrains_step(ntrls,a,b,T)
% INPUT:
% ntrls - number of trials
% a - initial firing rate (spikes/sec)
% b - final firing rate (spikes/sec)
% T - time period [seconds]
% OUTPUT
% s - spike train (ntrls*time ; 1ms bins) matrix containing zeros(no spike) and
%     ones(one spike) or twos(two spikes)
% n - step time [ms]

dt = 0.001;
lat_res = 10; % step latency resolution in ms

n=randi(T/(dt*lat_res),1,ntrls)*lat_res;
s=zeros(ntrls,T/dt);
for tr=1:ntrls
    for i= 1:n(tr)
        r=a;
        temp=rand;
        s(tr,i)=temp<=r*dt;
    end
    for i=n(tr)+1:T/dt
        r=b;
        temp=rand;
        s(tr,i)=temp<=r*dt;
    end
end
end