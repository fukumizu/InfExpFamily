%==========================================================================
%  run_score_bat_corr_new()
% 
%  run experiments score matching estimation with kernel exponential family
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

function run_score_bat_corr_new(n,d,type)


close all;

VERBOSE=0;
SAVEFILE=0;

nloop=50;
warning('off');

tau = 10.0; % standard deviation of base normal distribution


Datatype=type;

switch(Datatype)
    case 'Gauss'
        % True density: normal distriubtion
        sig_x=1.0;      % standard deviation 
        mean_x =0.0.*ones(1,d);
        param={mean_x, sig_x};
    case 'Gmix'
        K=2;    % number of components
        coefs=ones(K,1)./K; % coefficients for components (sum =1)
        means_x(1,:)=-4.*ones(1,d);  % means of components
        means_x(2,:)=4.*ones(1,d);
        sigs_x=1.0.*ones(K,1);      % sd of components (scalar matrix)
        param={coefs,means_x,sigs_x};
    otherwise
        error('Error: data type');
end
        
% Kernel parameters

% Kernel: exp(-||x-y||^2/(2*sx^2) + s (x^Ty + c)^2
% parameters 
c=0.5;
s=0.1;


if VERBOSE
    fprintf('Score Matching:\n');
    fprintf('kernel: c=%f, s=%f\n', c,s);
    fprintf('dim = %d, n = %d\n', d, n);
end

cor=zeros(nloop,1);
obj=zeros(nloop,1);
for loop=1:nloop
    % Initialization of random number generator 
    sd = RandStream('mt19937ar','Seed',loop);
    RandStream.setGlobalStream(sd);

    % Generating Data
    X = GenData(n,Datatype,param);

    [obj(loop), cor(loop), qu, pt]=Score_est_corr_new(X,c,s,tau,Datatype,param);
    
    if VERBOSE
        fprintf('%d\t%e\t%e\n', loop, obj(loop), cor(loop));
    end
    
end

if VERBOSE
    fprintf('\n');
    fprintf('dim = %d, n = %d\n', d, n);
    fprintf('Obj fun: mean = %.10f, std = %.10f\n', mean(obj), std(obj));    
    fprintf('Cor: mean = %f, std = %f\n', mean(cor), std(cor));
end


fprintf('%d\t%d\t%f\t%f\t%e\t%e\n',n,d, mean(obj), std(obj), mean(cor), std(cor));

if SAVEFILE
    fname=sprintf('result_n%dd%d%s',n,d,type);
    save(fname, 'obj','cor','qu','pt');
end




    
