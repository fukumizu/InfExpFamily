% test batch
%  used for the experiments in the following paper
%  "Density Estimation in Infinite Dimensional Exponential Families" 
%   by Bharath Sriperumbudur, Kenji Fukumizu, Arthur Gretton, Aapo
%   Hyvarinen and Revant Kumar, Journal of Machine Learning Research 2017 
%--------------------------------------------------------------------------

% Gauss score

fprintf('\nScore Matching Gauss (sample size)\n')
parfor i=1:5
    run_score_bat_corr_new(i*100,7,'Gauss');
end

fprintf('\nScore Matching Gauss (dimension)\n')
parfor i=1:10
    run_score_bat_corr_new(500,i*2,'Gauss');
end


fprintf('\nKDE Gauss (sample size)\n');

% Gauss KDE
for i=1:5
    run_kde_bat_corr(i*100,7,'Gauss')
end

fprintf('\nKDE Gauss (sample size)\n');
for i=1:10
    run_kde_bat_corr(500,2*i,'Gauss');
end


% Gmix Score

fprintf('\nScore Matching Gmix (sample size)\n')

parfor i=1:5
    run_score_bat_corr_new(100*i,7,'Gmix')
end

fprintf('\nScore Matching Gmix (dimension)\n')
parfor i=1:10
    run_score_bat_corr_new(300,2*i,'Gmix')
end

% Gmix KDE
fprintf('\nKDE Gmix (sample size)\n');
for i=1:5
    run_kde_bat_corr(100*i, 7, 'Gmix')
end
fprintf('\nKDE Gmix (dimension)\n');
for i=1:10
    run_kde_bat_corr(300,2*i,'Gmix');
end


