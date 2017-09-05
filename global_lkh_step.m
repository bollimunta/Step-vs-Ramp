function [GLK_step]=global_lkh_step(s,T,A,B)
% global_lkh_step(s,T,A,B) estimates global likelihood given that firing
% rate is stepping
%
% INPUT:
% s - spike train (ntrls*time ; 1ms bins) matrix containing zeros(no spike) and
%     ones(one spike) or twos(two spikes)
% T - time period in seconds (ex. 0.3 s of presaccadic activity)
% A - initial firing rate (spk/s) estimate
% B - final firing rate (spk/s) estimate
% OUPUT:
% GLK_step - global likelihood
% revised March 2017 by Anil Bollimunta

%% 
% compute log likelihood for each combination of
% parameters given that firing rate is stepping
[lk_step]=glh_step(s,T,A,B);

%% 
% convert loglikelihood to likelihood and sum over the parameter space
ntrls=size(lk_step,1);
Na=size(lk_step,2);
Nb=size(lk_step,3);
nbins=size(lk_step,4);

G=zeros(1,1);
for n1=1:nbins
    tr=1;
    x1=lk_step(tr,:,:,n1);
    for n2=1:nbins
        tr=2;
        x2=x1+lk_step(tr,:,:,n2);
        for n3=1:nbins
            tr=3;
            x3=x2+lk_step(tr,:,:,n3);
            for n4=1:nbins
                tr=4;
                x4=x3+lk_step(tr,:,:,n4);
                    G = G + sum(sum(exp(x4))); % if computing over 4 trials uncomment this line and remove the 5th for-loop
%                 for n5=1:nbins
%                     tr=5;
%                     x5=x4+lk_step(tr,:,:,n5);
%                     G = G + sum(sum(exp(x5)));
%                 end
            end
        end
    end
end
temp=nbins^ntrls;
GLK_step =G/(Na*Nb*temp);

end

%%
function [lk_step]=glh_step(s,T,A,B)
% glh_step(s,T,A,B) computes log likelihood for each combination of
% parameters given that firing rate is stepping

% set firing rate(spk/s) parameter range with 2 Hz resolution
a=A-50:2:A+50;
x=a>0;a=a(x);
b=B-50:2:B+50;
y=b>0;b=b(y);

% step latency resolution in ms
lat_res=10; 

dt=0.001; % millisecond 
n=1:lat_res:T/dt;
ntrls=size(s,1);
nbins=size(n,2);
lk_step=zeros(ntrls,size(a,2),size(b,2),nbins);
for tr=1:ntrls
    spk=s(tr,:);
    for i=1:size(a,2)
        for j=1:size(b,2)
            for k=1:nbins
                lk_step(tr,i,j,k)=loglkhd_step_1trl(spk,T,a(i),b(j),k*lat_res);
            end
        end
    end
end
% mx_lkd=max(max(max(lk_step)));
end
%%
function Q=loglkhd_step_1trl(s,T,a,b,lat) 
% loglkhd_step_1trl(s,T,a,b,lat) estimates log likelihood of a single trial given that firing
% rate is stepping

dt=0.001;
[x,y]=find(s==1);[x2 y2]=find(s==2);

rate(1:lat)=a;
rate(lat+1:T/dt)=b;

Q=-sum(rate)*dt + sum(log(rate(y))) + 2*sum(log(rate(y2)));
end
