%%
clear all
close all

%%
r=.1; sigma=.4; K=1;
S_range=0:.05*K:2*K;
tau_range=0:.05:3;
[S,tau]=meshgrid(S_range,tau_range);
[c,p,deltac,deltap,gamma,thetac,thetap,vega,volga] = optionCalc(S,tau,r,sigma,K);

figure()
hold on
mesh(S,tau,c)
mesh(S,tau,p);set(gca,'FontSize',14)
hold off
xlabel('S/K')
ylabel('\tau')
zlabel('c, p')

% figure()
% hold on
% mesh(S,tau,deltac)
% mesh(S,tau,deltap)
% hold off
% xlabel('S/K')
% ylabel('\tau')
% zlabel('\Delta_c, \Delta_p')
% 
% figure()
% mesh(S,tau,gamma)
% xlabel('S/K')
% ylabel('\tau')
% zlabel('\Gamma')
% 
% figure()
% hold on
% mesh(S,tau,thetac)
% mesh(S,tau,thetap)
% hold off
% xlabel('S/K')
% ylabel('\tau')
% zlabel('\Theta_c, \Theta_p')

figure()
mesh(S,tau,vega);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\upsilon')

figure()
mesh(S,tau,volga);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\partial \upsilon / \partial \sigma')

cp_parity=c+K*exp(-r*tau)-p-S;
disp(['check for call-put parity: ' num2str(sum(sum(cp_parity.^2)))])
cp_delta_diff=deltac-deltap;

%%
figure()
subplot(2,2,1)
mesh(S,tau,c);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('c')
subplot(2,2,3)
mesh(S,tau,deltac);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\Delta_c')
subplot(2,2,2)
mesh(S,tau,p);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('p')
subplot(2,2,4)
mesh(S,tau,deltap);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\Delta_p')

figure()
subplot(2,2,1)
mesh(S,tau,c);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('c')
subplot(2,2,3)
mesh(S,tau,gamma);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\Gamma')
subplot(2,2,2)
mesh(S,tau,p);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('p')
subplot(2,2,4)
mesh(S,tau,gamma);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\Gamma')

figure()
subplot(2,2,1)
mesh(S,tau,c);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('c')
subplot(2,2,3)
mesh(S,tau,thetac);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\Theta_c')
subplot(2,2,2)
mesh(S,tau,p);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('p')
subplot(2,2,4)
mesh(S,tau,thetap);set(gca,'FontSize',14)
xlabel('S/K')
ylabel('\tau')
zlabel('\Theta_p')