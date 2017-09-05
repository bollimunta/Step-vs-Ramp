function [GLK_ramp] = global_lkh_ramp(s,T,A,B)
% global_lkh_ramp(s,T,A,B) estimates global likelihood given that firing
% rate is ramping
%
% INPUT:
% s - spike train (ntrls*time ; 1ms bins) matrix containing zeros(no spike) and
%     ones(one spike) or twos(two spikes)
% T - time period in seconds (ex. 0.3 s of presaccadic activity)
% A - initial firing rate (spk/s) estimate
% B - final firing rate (spk/s) estimate
% OUPUT:
% GLK_ramp - global likelihood
% revised March 2017 by Anil Bollimunta

%%
% set firing rate(spk/s) parameter range with 2 Hz resolution
a=A-50:2:A+50; 
x=a>0;a=a(x);
b=B-50:2:B+50;
y=b>0;b=b(y);

%% 
Na=size(a,2);
Nb=size(b,2);
lk_ramp=zeros(size(a,2),size(b,2));
for i=1:size(a,2)
    for j=1:size(b,2)
        lk_ramp(i,j)=loglkhd_ramp(s,a(i),b(j),T);
    end
end

G=sum(sum(exp(lk_ramp))); % convert loglikelihood to likelihood and sum over the parameter space
GLK_ramp=G/(Na*Nb);
% mx_lkd=max(max(lk_ramp)); % max likelihood

end

function Q=loglkhd_ramp(s,a,b,T) 
% loglkhd_ramp(s,a,b,T) estimates log likelihood given that firing
% rate is ramping  
ntrls=size(s,1);
dt=0.001; % millisecond bin size
bin=1:T/dt;

[x,y]=find(s==1);[x2, y2]=find(s==2);
rate=a + (b-a)*bin/size(bin,2);

% Q = -ntrls*(a*T+ b*T*T/2) + sum(log(rate(y)));
  
Q = -ntrls*(a*T + (b-a)*T/2) + sum(log(rate(y))) + 2*sum(log(rate(y2)));
end



