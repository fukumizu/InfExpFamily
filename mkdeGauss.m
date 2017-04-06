function f=mkdeGauss(trainX, testX, sx)

[ntr,d]=size(trainX);
[ntest,d]=size(testX);

UX=testX*trainX';
UU=sum(testX.*testX,2);
XX=sum(trainX.*trainX,2);
Dux=-2.*UX+repmat(UU,1,ntr)+repmat(XX',ntest,1);
f=mean(exp(-Dux./(2*sx*sx)),2);
f=f./(sqrt(2*pi)*sx)^d;