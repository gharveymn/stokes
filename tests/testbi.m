n=10; n1=n-1; n12=n1^2; h=1/n;
A=n^2*(2*diag(ones(n1,1))-diag(ones(n-2,1),1)-diag(ones(n-2,1),-1)); %1d laplacian
I=eye(n1); C=zeros(n1,n1); C(1,1)=2*n^4; C(n1,n1)=2*n^4; Z = zeros(n1,n1);
B = kron(A^2,I)+kron(I,A^2)+kron(A,A)+kron(C,I)+kron(I,C);
f=ones(n12,1);
lap2 = laplacian2(n1,n1,h);
M = [-lap2 (kron(I,I)+lap2)
	kron(Z,Z)    -lap2];
u1 = M\[f;f];
u1 = u1(1:n12);
u2 = B\f;
norm(u1-u2,inf)
U1 = reshape(u1,[n1,n1])';
U2 = reshape(u2,[n1,n1])';