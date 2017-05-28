% solve using discrete separation of variables   AUA+AUC+CUA=F, A,C sym
% Generalized eigval problem CP=APD, D diagonal, P orthogonal
% I=P^T A P, D=P^T C P. Then V=P^T U P. Then V+VD+DV  = P^TFP = G
% Then V+VD+DV= G or V_{ij}(1+d_i+d_j)= G_{ij}
n=10; n1=n-1; n12=n1^2; h=1/n;
A=n^2*(2*diag(ones(n1,1))-diag(ones(n-2,1),1)-diag(ones(n-2,1),-1)); %1d laplacian
I=eye(n1); C=zeros(n1,n1); C(1,1)=2*n^4; C(n1,n1)=2*n^4;
AA=kron(A,A)+kron(C,A)+kron(A,C);
[P,D]=eig(C,A);
f=ones(n12,1);
%G=(P\(A\reshape(f,n1,n1)))*P;
G=P'*reshape(f,n1,n1)*P;
V=G./(1+repmat(diag(D),1,n1)+repmat(diag(D)',n1,1));
%y=reshape(P*((Z/P)/A),n1^2,1);
%y=reshape(P'\(Z/P),n12,1);
U = P*V*P';
y=reshape(U,n12,1);
norm(y-AA\f,inf)