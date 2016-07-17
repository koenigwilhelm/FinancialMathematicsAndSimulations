%%
close all
clear all

%%
S0=40;
K=40;
T=2;
r=.06;
sigma=.4;
npath=10000;
ntime=1000;
dt=T/ntime;
tt=0:dt:T;

paths=zeros(npath,ntime+1)+S0;
coef=r-.5*sigma^2;
for ix=1:ntime
    zrnd=randn(npath,1);
    paths(:,ix+1)=paths(:,ix).*exp(coef*dt+sqrt(dt)*sigma*zrnd);
end

% figure()
% plot(tt'*ones(1,npath),paths')

V=zeros(npath,ntime+1);
V(:,end)=max(K-paths(:,end),0);
for ix=ntime:-1:1
    excpayoff=max(K-paths(:,ix),0);
    V(:,ix)=V(:,ix+1)*exp(-r*dt);
    if (ix>1)
        cregrs=polyfit(paths(excpayoff>0,ix),V(excpayoff>0,ix),2);
        indsel=(excpayoff>cregrs(1)*paths(:,ix).^2+cregrs(2)*paths(:,ix)+cregrs(3))&(excpayoff>0);
    else
        indsel=excpayoff>V(:,ix);
    end
    V(indsel,ix)=excpayoff(indsel);
end
disp(['put option price= ' num2str(mean(V(:,2))*exp(-r*dt))])