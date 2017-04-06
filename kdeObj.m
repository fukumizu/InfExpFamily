function fobj=kdeObj(X,U,sx,Datatype,Param)

[n d]=size(X);
neval=length(U(:,1));
sx2=sx*sx;

%
% -- computation of d logp(U)/dx
UX=U*X';
UU=sum(U.*U,2);
XX=sum(X.*X,2);
Dux=-2.*UX+repmat(UU,1,n)+repmat(XX',neval,1);
Kux=exp(-Dux./(2*sx2))./(sqrt(2*pi)*sx)^d;

p=mean(Kux,2);
H=-U.*repmat(p,1,d)+(Kux*X)./n;
H=H./repmat(p,1,d)./sx2;

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

fobj=sum(sum((H-H0).^2,1),2)/(2*neval);

    
    
    


