function [A,B]=Nodal_param(Z,Y,long)

[nc,nc,Ns]=size(Z);

for i=1:Ns
    [M,D]=eig(squeeze(Z(:,:,i))*squeeze(Y(:,:,i)));
    Mi=inv(M);
    Landa=sqrt(D);
    Yo=inv(squeeze(Z(:,:,i)))*M*Landa*Mi;
    Zo=inv(Yo);
    A(:,:,i)=Yo*M*diag(coth(diag(Landa)*long))*Mi;                                                 % Based on trigonometric Functions
    B(:,:,i)=-Yo*M*diag(csch(diag(Landa)*long))*Mi;                                                 % Based on trigonometric Functions
%     A(:,:,i)=Yo*(M*diag((ones(nc,1)+exp(-2*diag(D)))./(ones(nc,1)-exp(-2*diag(D))))*(Mi));  % Based on exponential Functions
%     B(:,:,i)=Yo*(M*diag((2*exp(-diag(D)))./(ones(nc,1)-exp(-2*diag(D))))*Mi);               % Based on exponentials Functions
end