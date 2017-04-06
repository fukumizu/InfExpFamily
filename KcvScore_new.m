%==========================================================================
%  KcvScore_new()
% 
%  Cross-validation for score mathcing with kernel expoenential family 
%
%  Author: Kenji Fukumizu
%  Affiliation: The Institute of Statistical Mathematics, ROIS
%  Date: April 6, 2017
%  Version: 1.00
%  Copyright: Kenji Fukumizu
%------------------------------------------------------------------------
%  
%  This matlab code calls the main routine for computing the score matching
%  (unnormalized) density estimator with infinte dimensional kernel
%  exponential family.  See the following paper for the details.
%
%  "Density Estimation in Infinite Dimensional Exponential Families" 
%   by Bharath Sriperumbudur, Kenji Fukumizu, Arthur Gretton, Aapo
%   Hyvarinen and Revant Kumar, Journal of Machine Learning Research 2017 
%
%==========================================================================

function [sx, lambda, CVerror]=KcvScore_new(X,sx_tbl,lam_tbl,c,s,tau,K)

VERBOSE = 0;

[n, d]=size(X);
tau2=tau*tau;

rX=X(randperm(n),:);

CVerror=zeros(length(sx_tbl),length(lam_tbl));
nk=floor(n/K);
for h=1:length(sx_tbl)
    sx=sx_tbl(h);
    sx2=sx*sx;
    for j=1:length(lam_tbl)
        lambda=lam_tbl(j);
        Kerrors=zeros(K,1);
        for k=1:K
            endp=min(k*nk,n);
            test=(k-1)*nk+1:endp;
            train=1:n;
            train(test)=[];
            trainX=rX(train',:);
            testX=rX(test,:);
            
            ntr=length(train);
            ntest=length(test);
            
            [alpha, beta] = ScoreMatch_new(trainX,sx,c,s,tau,lambda);
            beta=reshape(beta,ntr,d);
            
            % Gram matrix for testX and trainX 
            UX=testX*trainX';
            UU=sum(testX.*testX,2);
            XX=sum(trainX.*trainX,2);
            Dux=-2.*UX+repmat(UU,1,ntr)+repmat(XX',ntest,1);
            Kux=exp(-Dux./(2*sx2));
            
            Ub=testX*beta';
            Xb=sum(trainX.*beta,2);
            
            % d^2 f(X_t) / dx^2
            %Gauss
            d2f=mean(mean(Kux,1),2)*d*(d+2)/sx2/sx2;
            d2f=d2f-mean(mean(Dux.*Kux,1),2)*(2*d+4)/(sx2^3);
            d2f=d2f+mean(mean(Dux.*Dux.*Kux,1),2)/(sx2^4);
            Mua=(UX-repmat(XX',ntest,1));
            M=Mua.*((d+2)/sx2/sx2/tau2);
            M=M-Dux.*Mua./(tau2*sx2^3);
            d2f=d2f+mean(mean(M.*Kux,1),2);
            clear M;
            d2f=d2f*alpha;
            M2=Ub-repmat(Xb',ntest,1);
            d2f=d2f-mean(sum(M2.*Kux,2),1)*(d+2)/sx2/sx2;
            d2f=d2f+mean(sum(Dux.*M2.*Kux,2),1)/sx2^3;
            %poly
            d2fp=alpha*4*(d-mean(XX,1)/tau2) + 4*sum(Xb,1);
            %total
            d2ft=d2f+s*d2fp;
            
     
            % df/dx * dlogq/dx
            %Gauss
            Muu=(repmat(UU,1,ntr)-UX);
            R=-(d+2).*Muu./(tau2*sx2*sx2);
            R=R+Dux.*Muu./(tau2*sx2^3);
            R=R+UX./(tau2*tau2*sx2);
            R=R-Muu.*Mua./((sx2^2)*(tau2^2));
            dfq=mean(mean(R.*Kux,1),2)*alpha;
            dfq=dfq-mean(sum(Ub.*Kux,2),1)/tau2/sx2;
            dfq=dfq+mean(sum(Muu.*(Ub-repmat(Xb',ntest,1)).*Kux,2),1)/tau2/sx2^2;
            %poly
            dfqp=-4*mean(mean(UU,1),2)/tau2 + mean(mean(UX.*UX,1),2)*(4/tau2/tau2);
            dfqp=dfqp+mean(mean(UX,1),2).*(2*c/tau2/tau2);
            dfqp=dfqp*alpha;
            dfqp=dfqp - 4*mean(sum(Ub.*UX,2),1)/tau2;
            dfqp=dfqp - 2*c*mean(sum(Ub,2),1)/tau2;
            %total
            dfqt = dfq+s*dfqp;
            
            % df/dx
            %Gauss
            H=(d+2).*testX.*repmat(mean(Kux,2),1,d)./(sx2^2);
            H=H-Kux*trainX.*((d+2)/ntr/sx2/sx2);
            H=H-repmat(mean(Dux.*Kux,2),1,d).*testX./sx2^3;
            H=H+(Dux.*Kux)*trainX./ntr/sx2^3;
            H=H-Kux*trainX./ntr/sx2/tau2;
            H=H+repmat(mean(Mua.*Kux,2),1,d).*testX./sx2/sx2/tau2;
            H=H-(Mua.*Kux)*trainX./ntr/sx2/sx2/tau2;    % bug fixed: sept 24, Kenji
            H = H.*alpha;
            M3=M2.*Kux;
            H=H + Kux*beta./sx2 - repmat(sum(M3,2),1,d).*testX./sx2^2 ...
                + M3*trainX./sx2^2;
            %Poly
            Hp = 4.*testX - UX*trainX.*(4/ntr/tau2);
            Hp = Hp - 2.*c.*repmat(mean(trainX,1),ntest,1)./tau2;
            Hp=Hp*alpha;
            Hp = Hp + 2.*Ub*trainX + 2.*UX*beta + (2*c).*repmat(sum(beta,1),ntest,1);  
            %total
            Ht = H+s.*Hp;
            
            Kerrors(k)=sum(mean(Ht.*Ht,1),2)/2 + d2ft + dfqt;
        end
        CVerror(h,j)=mean(Kerrors);
    end
end

[minlam, isx]=min(CVerror);
[minval, lamidx]=min(minlam);

lambda=lam_tbl(lamidx);
sx=sx_tbl(isx(lamidx));

if VERBOSE
    fprintf('\nmin=%f,  lam = %f (%d/%d),  sx = %f (%d/%d)\n', ...
        minval,lambda, lamidx,length(lam_tbl), sx, isx(lamidx),length(sx_tbl));
end



            
            
            
            
            

            
            
            
            
            
        