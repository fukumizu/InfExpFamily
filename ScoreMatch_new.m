%==========================================================================
%  ScoreMatch_new()
% 
%  Main algorithm of score matching estimation for kernel exponential family
%
%  Author: Kenji Fukumizu
%  Affiliation: The Institute of Statistical Mathematics, ROIS
%  Date: April 6, 2017
%  Version: 1.00
%  Copyright: Kenji Fukumizu
%------------------------------------------------------------------------
%  
%  This matlab code is the main routine for computing the score matching
%  (unnormalized) density estimator with infinte dimensional kernel
%  exponential family.  See the following paper for the details.
%
%  "Density Estimation in Infinite Dimensional Exponential Families" 
%   by Bharath Sriperumbudur, Kenji Fukumizu, Arthur Gretton, Aapo
%   Hyvarinen and Revant Kumar, Journal of Machine Learning Research 2017 
%
%==========================================================================


function [alpha, beta] = ScoreMatch_new(X,sx,c,s,tau,lambda)

[n,d]=size(X);
tau2=tau*tau;
sx2=sx*sx;
ab=X*X';
aa=diag(ab);
D=repmat(aa,1,n);
Diff2=(D + D' - 2*ab);
K=exp(-Diff2./(2*sx2));    % Gauss Gram matrix

% G^{ab}_{ij} for Gaussian kernel
G=kron(eye(d),K./sx2);
G=reshape(G,n,d,n,d);
G=permute(G,[1 3 2 4]);     % G: nxnxdxd 

M=permute(repmat(X,[1,1,n]),[1 3 2]);
M=M-permute(M,[2,1,3]);     % M = (X^a_i - X^b_i)

G2=repmat(K./(sx2*sx2),[1 1 d d]);
B=repmat(M,[1,1,1,d]);
G2 = G2.*B;
G2 = G2.*permute(B,[1 2 4 3]);
G = G-G2;       
clear G2;

% G for polynomial kernel
Gp=kron(eye(d),2.*(ab+c));
Gp=reshape(Gp,n,d,n,d);
Gp=permute(Gp,[1 3 2 4]);     % G: nxnxdxd 
Ma=permute(repmat(X,[1 1 n d]),[1 3 4 2]);
Mb=permute(Ma,[2 1 4 3]);
Gp=Gp+2.*Ma.*Mb;

% G for Gauss K + s.polyn
Gt = G+s.*Gp;
% \sum_cm G_im^ac G_mj^cb
if d==1
    GG=Gt*Gt';
else
    A = reshape(permute(Gt,[1 3 2 4]),n*d,n*d);
    GG = A*A';
    clear A
end

% h vector for Gaussian kernel
h1 = -K*X;
h1=h1./(sx2*tau2);
h12=(repmat(aa,1,n)-ab).*K;
h12=h12./(sx2*sx2*tau2);
h12 = h12'*X - repmat(sum(h12,1)',1,d).*X;
h1=(h1+h12)./n;
h21=-K*X + repmat(sum(K,2),1,d).*X;
h21=h21.*((d+2)/(sx2*sx2));
h22=K.*Diff2./(sx2^3);
h2 = (h21 + h22*X - repmat(sum(h22,2),1,d).*X)./n;
h=h1+h2;

% h vector for polyn kernel
hp=-(ab*X).*(4/n/tau2) - repmat(mean(X,1),n,1).*(2*c/tau2) + 4.*X;

% h vector for Gauss + s.*polyn
ht=h+s.*hp;


% Solution (\alpha,\beta)
alpha = -1/lambda;  % In the new method, alpha = -1/\lambda;
Gt=reshape(permute(Gt,[1 3 2 4]),n*d,n*d);
beta = (Gt+n*lambda*eye(n*d))\(reshape(ht,n*d,1)/lambda);

