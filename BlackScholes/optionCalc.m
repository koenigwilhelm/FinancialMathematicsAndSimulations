function [c,p,deltac,deltap,gamma,thetac,thetap,vega,volga] = optionCalc(S,tau,r,sigma,K)
    % building blocks
    lnKdS=log(K./S);
    d1=((r+.5*sigma^2)*tau-lnKdS)./(sigma*sqrt(tau));
    d2=((r-.5*sigma^2)*tau-lnKdS)./(sigma*sqrt(tau));
    d1(isnan(d1))=0;
    d2(isnan(d2))=0;
    % call/put option price
    c=S.*normcdf(d1)-K*exp(-r*tau).*normcdf(d2);
    p=-S.*normcdf(-d1)+K*exp(-r*tau).*normcdf(-d2);
    % call/put deltas
    deltac=normcdf(d1);
    deltap=-normcdf(-d1);
    % call/put gamma
    gamma=exp(-d1.^2/2)./(sigma*S.*sqrt(2*pi*tau));
    % call/put thetas
    thetac=-S.*exp(-d1.^2/2)*sigma./(2*sqrt(2*pi*tau))-r*K*exp(-r*tau).*normcdf(d2);
    thetap=-S.*exp(-d1.^2/2)*sigma./(2*sqrt(2*pi*tau))+r*K*exp(-r*tau).*normcdf(-d2);
    % call/put vega and volga
    vega=S.*sqrt(tau/(2*pi)).*exp(-.5*d1.^2);
    volga=vega.*d1.*d2/sigma;
end