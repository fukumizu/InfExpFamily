% run experiments score matching estimation with kernel exponential family

function [fobj corr pu pt]=kde_est_corr(X,Datatype,param)

K=5;  % K-fold CV
[n,d]=size(X);

% K-fold CV to choose bandwidth of kernel
bw0=MedianDist(X);
bw=[0.02 0.04 0.06 0.08 0.1 0.2 0.4 0.6 0.8 1.0].*bw0;

rX=X(randperm(n),:);
CVloglik=zeros(length(bw),1);
nk=floor(n/K);

for h=1:length(bw)
    loglik=zeros(K,1);
    for k=1:K
        endp=min(k*nk,n);
        test=(k-1)*nk+1:endp;
        train=1:n;
        train(test)=[];
        trainX=rX(train',:);
        testX=rX(test,:);

        %loglik(k)=sum(log( ksdensity(trainX,testX,'width',bw(h))));
        loglik(k)=sum(log(mkdeGauss(trainX,testX,bw(h))));
    end
    CVloglik(h)=mean(loglik);
end
[minv, idx]=max(CVloglik);
optbw=bw(idx);
%fprintf('cv idx = %d\n', idx);
%fprintf('CV: width = %f (%d/%d)\n', optbw, idx, length(bw));

% Evaluation of objective function with new data points. 
%testX = GenData(5*n,Datatype,param);    % test data
%z=ksdensity(trainX,testX,'width',optbw);
%fprintf('Testing objective funciton = %f\n\n', z);


ntest=10000;  
U = GenData(ntest,Datatype,param);

pu=mkdeGauss(trainX,U,optbw);
pt = eval_true_density(U,Datatype,param);
corr=pu'*pt/sqrt(pu'*pu)/sqrt(pt'*pt);
fobj=kdeObj(U,trainX,optbw,Datatype,param);
%fprintf('KL div = %f   MSE = %.15f\n',div, diff2);



    




    
