%==========================================================================
%  run_score_bat_corr_new()
% 
%  run experiments for kernel density estimation
%
%  Author: Kenji Fukumizu
%  Affiliation: The Institute of Statistical Mathematics, ROIS
%  Date: April 6, 2017
%  Version: 1.00
%  Copyright: Kenji Fukumizu
%------------------------------------------------------------------------
%  
%  This matlab code calls the main routine for computing the standard
%  kernel density estimation for comparison.  See the following paper for
%  the details. 
%
%  "Density Estimation in Infinite Dimensional Exponential Families" 
%   by Bharath Sriperumbudur, Kenji Fukumizu, Arthur Gretton, Aapo
%   Hyvarinen and Revant Kumar, Journal of Machine Learning Research 2017 
%
%==========================================================================

function run_kde_bat_corr(n,d,type)

VERBOSE=0;
SAVEFILE=0;

nloop=50;


% Datatype='Gauss';   % Gaussian
% Datatype='Gmix';  % Gaussian mixture 
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


if VERBOSE
    fprintf('KDE:\n');
    fprintf('kernel: gauss\n');
    fprintf('dim = %d, n = %d\n', d, n);
end


corkde=zeros(nloop,1);
scorekde=zeros(nloop,1);
for loop=1:nloop
    % Initialization of random number generator 
    sd = RandStream('mt19937ar','Seed',loop);
    RandStream.setGlobalStream(sd);

    % Generating Data
    X = GenData(n,Datatype,param);

    [fobj, corr, pu, pt]=kde_est_corr(X,Datatype,param);

    
    if VERBOSE
        fprintf('%d\t%f\t%f\n', loop, fobj, corr);
    end
    
    corkde(loop)=corr;
    scorekde(loop)=fobj;
   
end

if VERBOSE
    fprintf('\nKDE: data = %s\n',Datatype)
    fprintf('N = %d   dim = %d', n, d)
    fprintf('\n');
    fprintf('dim = %d, n = %d\n', d, n);
    fprintf('Obj fun: mean = %.10f, std = %.10f\n', mean(scorekde), std(scorekde));    
    fprintf('Cor: mean = %f, std = %f\n', mean(corkde), std(corkde));
end


fprintf('%d\t%d\t%f\t%f\t%f\t%f\n',...
    n,d, mean(scorekde), std(scorekde), mean(corkde), std(corkde));

if SAVEFILE
    fname=sprintf('result_kde_n%dd%d%s',n,d,Datatype);
    save(fname, 'scorekde','corkde','pu','pt');
end
    
