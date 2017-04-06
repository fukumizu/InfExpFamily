# InfExpFamily
Matlab Codes for score matcning density estimation with infinite dimensional kernel exponential family 

The codes were used for generating the experiments in the following paper:
  Density Estimation in Infinite Dimensional Exponential Families
  by Bharath Sriperumbudur, Kenji Fukumizu, Arthur Gretton, Aapo Hyvarinen and Revant Kumar
  Journal of Machine Learning Research 2017

The results in the paper are given by test_score_bat_new.m (it takes much time).
To see a result of score matching estimation, use run_score_bat_corr_new.
  
  run_score_bat_corr_new(n,d,type)
    n: sample size
    d: dimensionality
    type: "Gauss" or "Gmix"  (Gmix is Gaussian mixture of two components)
    
Kernel density estimation can be done with run_kde_bat_corr(n,d,type).



