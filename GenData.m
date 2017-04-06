function X = GenData(n,Datatype,param)

switch(Datatype)
    case 'Gauss'
        mx=param{1};
        sx=param{2};
        d=size(mx,2);
        X=repmat(mx,n,1)+sx.*randn(n,d);
    case 'Gmix'
        coefs=param{1};
        mx=param{2};
        sx=param{3};
        cc=cumsum(coefs);
        d=size(mx,2);
        K=size(coefs,1);
        cls=rand(n,K)<repmat(cc',n,1);
        cls=cls-[zeros(n,1) cls(:,1:K-1)];
        X=cls*mx + repmat((cls*sx),1,d).*randn(n,d);
    otherwise 
        error('Error: Data type');
end

        
        
        