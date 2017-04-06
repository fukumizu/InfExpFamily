function pt = eval_true_density(U,Datatype,param)

[n,d]=size(U);
switch(Datatype)
    case 'Gauss'
        mean_x=param{1};
        sig_x=param{2};
        pt =exp( -sum((U-repmat(mean_x,n,1)).^2,2)./(2*sig_x*sig_x) )./((sqrt(2*pi)*sig_x)^d);
    case 'Gmix'
        coefs=param{1};
        mx=param{2};
        sx=param{3};
        d=size(mx,2);
        K=size(coefs);
        pt=zeros(n,1);
        for i=1:K
            pt=pt+coefs(i).*exp(-sum((U-repmat(mx(i,:),n,1)).^2,2)./(2*sx(i)*sx(i)))./(sqrt(2*pi)*sx(i))^d;
        end
    otherwise 
        error('Error: Data type');
end

        
        

