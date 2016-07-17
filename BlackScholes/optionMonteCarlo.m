function [c,p] = optionMonteCarlo(S,tau,r,sigma,K,numSim,flagImpSamp)
    [m,n]=size(S);
    c=zeros(m,n);
    p=c;
    nmr=normrnd(0,1,numSim,1);
    for ix=1:n
        scrateImpSamp=(log(K/S(1,ix))-.5*sigma^2)*tau(:,ix)*ones(1,numSim)+sigma*sqrt(tau(:,ix))*nmr';
        scrate=(r-.5*sigma^2)*tau(:,ix)*ones(1,numSim)+sigma*sqrt(tau(:,ix))*nmr';
        Sfin=S(:,ix)*ones(1,numSim).*exp(scrate);
        payoffC=max(Sfin-K,0);
        payoffP=max(K-Sfin,0);
        if (flagImpSamp==1)
            if (S(1,ix)<K)
                payoffC=max(S(:,ix)*ones(1,numSim).*exp(scrateImpSamp)-K,0);
                payoffC=payoffC.*normpdf(scrateImpSamp,(r-.5*sigma^2)*tau(:,ix)*ones(1,numSim),sigma*sqrt(tau(:,ix))*ones(1,numSim))...
                    ./normpdf(scrateImpSamp,(log(K/S(1,ix))-.5*sigma^2)*tau(:,ix)*ones(1,numSim),sigma*sqrt(tau(:,ix))*ones(1,numSim));
            elseif (S(1,ix)>K)
                payoffP=max(K-S(:,ix)*ones(1,numSim).*exp(scrateImpSamp),0);
                payoffP=payoffP.*normpdf(scrateImpSamp,(r-.5*sigma^2)*tau(:,ix)*ones(1,numSim),sigma*sqrt(tau(:,ix))*ones(1,numSim))...
                    ./normpdf(scrateImpSamp,(log(K/S(1,ix))-.5*sigma^2)*tau(:,ix)*ones(1,numSim),sigma*sqrt(tau(:,ix))*ones(1,numSim));
            end
            payoffC(isnan(payoffC))=0;
            payoffP(isnan(payoffP))=0;
        end
        c(:,ix)=mean(payoffC,2).*exp(-r*tau(:,ix));
        p(:,ix)=mean(payoffP,2).*exp(-r*tau(:,ix));
    end
end

