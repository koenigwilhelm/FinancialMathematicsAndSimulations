function [c,p] = optionMonteCarlo(S,tau,r,sigma,K,numSim)
    [m,n]=size(S);
    c=zeros(m,n);
    p=c;
    nmr=normrnd(0,1,numSim,1);
    for ix=1:n
        scrate=(r-.5*sigma^2)*tau(:,ix)*ones(1,numSim)+sigma*sqrt(tau(:,ix))*nmr';
        scrate=exp(scrate);
        Sfin=S(:,ix)*ones(1,numSim).*scrate;
        payoffC=Sfin-K;
        payoffP=K-Sfin;
        payoffC(payoffC<0)=0;
        payoffP(payoffP<0)=0;
        c(:,ix)=mean(payoffC,2).*exp(-r*tau(:,ix));
        p(:,ix)=mean(payoffP,2).*exp(-r*tau(:,ix));
    end
end

