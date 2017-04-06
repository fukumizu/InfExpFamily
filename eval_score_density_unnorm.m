% 
% function [pu, zu, fu]=eval_score_density_unnorm(X,sx,c,s,tau,alpha,beta,U)
%
% INPUT
%  X:   trainig data
%  sx:  sigma for Gaussian kernel exp( - || x-y ||^2 / (2*sx*sx) )
%  c:   constant in polynomial kernel (x^T y + c)^2
%  s:   coefficient of kernel sum: k(x,y) = Gaussan(x,y) + s.*Polyn(x,y)
%  tau: standard deviation of base normal distribution
%  alpha, beta:  estimated parameter
%  
% OUTPUT
%  pu: values of density function at testing points
%  zu: values of partition function
%  fu: exp(f(u))  (unnormalized density function estimated)


function qu=eval_score_density_unnorm(X,sx,c,s,tau,alpha,beta,U)

[n,d]=size(X);
sx2=sx*sx;
tau2=tau*tau;

n_test=length(U(:,1));

%
% -- computation of f(U)
aa=sum(X.*X,2);
ux=U*X';
uu=sum(U.*U,2);
Dux=repmat(uu,1,n)+repmat(aa',n_test,1)-2.*ux;
Kux=exp(-Dux./(2*sx2));

% \xi(U) for Gaussian
xi=-mean(Kux,2).*(d/sx2);
xi=xi + mean(Dux.*Kux,2)./(sx2*sx2);
xi=xi + mean((ux-repmat(aa',n_test,1)).*Kux,2)/(sx2*tau2);
% \xi(Y) for polyn
xip = 2.*sum(U.*U,2) - (mean(ux.*ux,2)+c.*mean(ux,2)).*(2/tau2);

xit=xi+s.*xip;

% \sum \beta d k(U,X)/d
bm=reshape(beta,n,d);
ub=(bm*U')';    % n_test x n
xb=sum(bm.*X,2);     % n x 1
fua=alpha*xit;
fub=  sum( (ub-repmat(xb',n_test,1)).*Kux, 2)./sx2 ...
    +s.*(2.*sum(ux.*ub,2) + 2.*c.*sum(ub,2));       % bug fixed Sept 25, fukumizu
fu=fua+fub;

qu=exp(fu);
