%%
close all
clear all

%%
S0=42;
K=40;
T=2;
r=.06;
sigma=.4;
npath=10000;
ntime=1000;
dt=T/ntime;
tt=0:dt:T;

paths=zeros(2*npath,ntime+1)+S0;
coef=r-.5*sigma^2;
for ix=1:ntime
    zrnd=randn(npath,1);
    % original paths
    paths(1:npath,ix+1)=paths(1:npath,ix).*exp(coef*dt+sqrt(dt)*sigma*zrnd);
    % antithetic variates
    paths(npath+1:end,ix+1)=paths(npath+1:end,ix).*exp(coef*dt-sqrt(dt)*sigma*zrnd);
end

figure()
plot(tt'*ones(1,2*npath),paths');set(gca,'FontSize',14)
xlabel('t')
ylabel('S(t)')
title({'Stock Price Paths';['T=' num2str(T) '; r=' num2str(r) '; \sigma=' num2str(sigma) '; S_0=' num2str(S0)]})

V=zeros(2*npath,ntime+1);
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
%     paths(indsel,ix+1:end)=paths(indsel,ix)*ones(1,ntime+1-ix);
end
disp(['American put option price= ' num2str(mean(V(:,2))*exp(-r*dt))])
[~,p,~,~,~,~,~,~,~]=optionCalc(S0,T,r,sigma,K);
disp(['European put option price= ' num2str(p)])
% figure()
% plot(tt'*ones(1,2*npath),paths')