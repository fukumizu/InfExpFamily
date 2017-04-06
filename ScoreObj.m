function fobj=ScoreObj(X,sx,c,s,tau,alpha,beta,Datatype,Param)

neval=10000;

[n,d]=size(X);
sx2=sx*sx;
tau2=tau*tau;

U=GenData(neval,Datatype,Param);    % neval data from the true probability

%
% -- computation of df(U)/dx
UX=U*X';
UU=sum(U.*U,2);
XX=sum(X.*X,2);
Dux=-2.*UX+repmat(UU,1,n)+repmat(XX',neval,1);
Kux=exp(-Dux./(2*sx2));
            
Ub=U*beta';
Xb=sum(X.*beta,2);

%Gauss
H=(d+2).*U.*repmat(mean(Kux,2),1,d)./(sx2^2);
H=H-Kux*X.*((d+2)/n/sx2/sx2);
H=H-repmat(mean(Dux.*Kux,2),1,d).*U./sx2^3;
H=H+(Dux.*Kux)*X./n/sx2^3;
H=H-Kux*X./n/sx2/tau2;
H=H+repmat(mean((UX-repmat(XX',neval,1)).*Kux,2),1,d).*U./sx2/sx2/tau2;
H=H-((UX-repmat(XX',neval,1)).*Kux)*X./n/sx2/sx2/tau2;    
H = H.*alpha;
M3=(Ub-repmat(Xb',neval,1)).*Kux;
H=H + Kux*beta./sx2 - repmat(sum(M3,2),1,d).*U./sx2^2 ...
                + M3*X./sx2^2;
clear M3;
%Poly
Hp = 4.*U - UX*X.*(4/n/tau2);
Hp = Hp - 2.*c.*repmat(mean(X,1),neval,1)./tau2;
Hp=Hp*alpha;
Hp = Hp + 2.*Ub*X + 2.*UX*beta + (2*c).*repmat(sum(beta,1),neval,1);  
%total
Hest = H+s.*Hp - U./tau2;

%
% -- computation of d log p_0(U) / dx

if strcmp(Datatype,'Gauss')
    means=Param{1};
    st=Param{2};
    H0 = -(U-repmat(means,neval,1))./(st*st);
elseif strcmp(Datatype,'Gmix')
    coef=Param{1};
    means=Param{2};
    st=Param{3};
    K=length(coef);
    H0=zeros(neval,d);
    p0=zeros(neval,1);
    for i=1:K
        M=U-repmat(means(i,:),neval,1);
        v=coef(i).*exp( -sum(M.*M,2)./(2*st(i)*st(i)))./(sqrt(2*pi)*st(i))^d;
        p0=p0+v;
        H0=H0-repmat(v,1,d).*M./(st(i)*st(i));
    end
    H0=H0./repmat(p0,1,d);
end

fobj=sum(sum((Hest-H0).^2,1),2)/(2*neval);

    
    
    


