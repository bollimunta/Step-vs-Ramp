function global_ratio_data(spk)
% global_ratio_LIP(spk) estimates global likelihood ratio of LIP spike
% trains and artificial(poisson) spike trains generated with similar
% parameters.
%
% INPUT:
% spk - spike train (ntrls*time ; 1ms bins) matrix containing zeros(no spike) and
%     ones(one spike) or twos(two spikes)
close all
%% initial and final firing rates
f = 1000; % sampling frequency (Hz)
bin_sz = 1; % bin size (ms)
a = mean(mean(spk(:,1:40)))*f/bin_sz;       % intial firing rate during 40 ms period (spk/s)
b = mean(mean(spk(:,end-40:end)))*f/bin_sz; % final firing rate during 40 ms period (spk/s)

%% run
T = size(spk,2)/f; % time period in seconds (ex. 0.3 s of presaccadic activity)
ntrls = size(spk,1); % number of trials
itns = 25; % number of iterations
ntr = 4; % number of trials over which global likelihhod ratio is computed during each iteration
w = waitbar(0,'iteration no');
for i= 1:itns
    %     display('computing Global LIP ratio')
    drw=randi(ntrls,1,ntr); % draw 'ntr' number of trials
    s=spk(drw,:);
    [GLK_ramp_LIP] = global_lkh_ramp(s,T,a,b);
    [GLK_step_LIP]=global_lkh_step(s,T,a,b);
    GLK_ratio_LIP(i)=GLK_ramp_LIP/GLK_step_LIP;
    
    clear s
    %     display('computing Global ramp ratio')
    [s]=poissontrains_ramp(ntr,a,b,T);
    [GLK_ramp_ramp] = global_lkh_ramp(s,T,a,b);
    [GLK_step_ramp]=global_lkh_step(s,T,a,b);
    GLK_ratio_ramp(i)=GLK_ramp_ramp/GLK_step_ramp;
    
    clear s
    %     display('computing Global step ratio')
    [s,~]=poissontrains_step(ntr,a,b,T);
    [GLK_ramp_step] = global_lkh_ramp(s,T,a,b);
    [GLK_step_step]=global_lkh_step(s,T,a,b);
    GLK_ratio_step(i)=GLK_ramp_step/GLK_step_step;
    waitbar(i/itns,w)
end
close(w)

L = log(GLK_ratio_LIP);
R = log(GLK_ratio_ramp);
S = log(GLK_ratio_step);

%% plot
h = figure(1);
set(h,'position',[-700 750 500 800],'color','w')
s1 = axes('position',[0.1 0.1 0.8 0.4],'Parent',h);
hold on
edges = linspace(min([L R S]),max([L R S]),25);
nL = histc(L,edges);
nR = histc(R,edges);
nS = histc(S,edges);
stairs(edges+(edges(2)-edges(1))/2,nL,'b','linewidth',2)
stairs(edges+(edges(2)-edges(1))/2,nR,'g','linewidth',2)
stairs(edges+(edges(2)-edges(1))/2,nS,'r','linewidth',2)
line([nanmedian(L) nanmedian(L)],[0 max([nL nR nS])+ 5],'linewidth',2,'linestyle','--','color','b')
line([nanmedian(R) nanmedian(R)],[0 max([nL nR nS])+ 5],'linewidth',2,'linestyle','--','color','g')
line([nanmedian(S) nanmedian(S)],[0 max([nL nR nS])+ 5],'linewidth',2,'linestyle','--','color','r')
grid on
set(get(s1,'Title'),'string','Step vs Ramp' ,'FontSize',12,'FontWeight','bold')
xlabel('Odds ratio log(ramp/step)','FontSize',12,'FontWeight','bold')
ylabel('no. of random draws','FontSize',12,'FontWeight','bold')
legend('LIP','ramp','step','Location','NorthWest')
legend('boxoff')

s2 = axes('position',[0.1 0.6 0.8 0.3],'Parent',h);
stairs(mean(spk)*f,'k')
grid on
set(get(s2,'Title'),'string','firing rate' ,'FontSize',12,'FontWeight','bold')
xlabel('time [ms]','FontSize',12,'FontWeight','bold')
ylabel('spikes/sec','FontSize',12,'FontWeight','bold')

%% save
temp = cd;
if ispc
    cd('Z:\System Backups\B-rig\Analysis\vprobe\cor')
else
    cd('/Users/bollimuntaak/Dropbox/poisson/basso')
end
set(gcf, 'PaperPositionMode', 'auto')   % Use screen size
set(gcf,'PaperUnits','normalized','PaperPosition',[0.2 0.1 0.6 0.8])
print(1, '-dpdf', 'GLK_LIP')
cd(temp)
save('log odds','GLK_ratio*','spk')

end

%% generate poisson spike trains with ramping fr 
function [s]=poissontrains_ramp(ntrls,a,b,T)

dt=0.001;
s=zeros(ntrls,T/dt);
for i= 1:T/dt
    r=a + (b-a)*i*dt/T;
    temp=rand(ntrls,1);
    s(:,i)=temp<=r*dt;
end
end

%% generate poisson spike trains with stepping fr 
function [s,n]=poissontrains_step(ntrls,a,b,T)

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
