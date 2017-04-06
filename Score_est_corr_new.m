%==========================================================================
%  Score_est_corr_new()
% 
%  Score matching estimation for kernel exponential family
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


function [obj, cor, qu, pt]=Score_est_corr_new(X,c,s,tau,Datatype,param)

CVK=5;  % K-fold CV
[n,d]=size(X);

DISP=0;


% K-fold CV to choose bandwidth of kernel and regularization coeff.
sx0=MedianDist(X);  % median of pairwise distances 
sx_tbl=[0.1 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6]*sx0;
lam_tbl=0.1/(n^(1/3));
[sx, lambda, CVerror]=KcvScore_new(X,sx_tbl,lam_tbl,c,s,tau,CVK);


[alpha, beta] = ScoreMatch_new(X,sx,c,s,tau,lambda);
beta=reshape(beta,n,d);


ntest=10000;
U = GenData(ntest,Datatype,param);

qu=eval_score_density_unnorm(X,sx,c,s,tau,alpha,beta,U);

pt = eval_true_density(U,Datatype,param);
cor=qu'*pt/sqrt(qu'*qu)/sqrt(pt'*pt);
obj=ScoreObj(X,sx,c,s,tau,alpha,beta,Datatype,param);


%figure;
if DISP
    plot(U,pu,'r.');
    hold on;
    plot(U,pt,'k.');
    hold off;
    drawnow;
end


    
    


    
