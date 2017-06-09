% finite difference Poisson solver on a square

nn=2.^(3:6); %[10,20,40,80,160]; 
for j=1:size(nn,2)
    n=nn(j); h=1/n; n1=n-1; N=n1*n1;
    x=h:h:1-h; [X,Y]=meshgrid(x,x);
    fn=inline('sin(pi*x).*sin(2*pi*y)','x','y');
    u=fn(X,Y); u=reshape(u,N,1); % exact soln
    fh=5*pi*pi*u; fh=reshape(fh,N,1); % RHS
    e=ones(N,1);
    A=spdiags([-e,-e,4*e,-e,-e],[-n1,-1,0,1,n1],N,N)/h^2;
    for kk=1:n-2 k=kk*n1; A(k,k+1)=0; A(k+1,k)=0; end;
    uh=A\fh; % computed soln by Gaussian elimination
    er(j)=norm(uh-u,inf);
end;
figure(1); plot(nn,er,'o-'); xlabel('n'); ylabel('error');
figure(2); surf(X,Y,reshape(uh,n1,n1)); xlabel('x'); ylabel('y'); zlabel('z');
return
% Solve L u = f where L is -Laplacian + random (sym) perturbation using PCG
clear X Y uh u
A=toeplitz([2 -1 zeros(1,n-3)])*n^2; % 1D Laplacian
[V,D]=eig(A); % spectral decomposition needed for tensor product preconditioner
A=kron(A,eye(n1))+kron(eye(n1),A)+toeplitz(randn(1,N)); % calculate L

[x,flag,relres,it,resv]=pcg(A,fh,[],[],@(f)tenpre(f,V,D,n1));
figure(3); semilogy([1:it]',[resv(1:it)],'-o'); xlabel('iteration #'); ylabel('error');

function u=tenpre(f,V,D,n) % tensor product preconditioner
G=V'*reshape(f,n,n)*V;
Z=G./(repmat(diag(D),1,n) + repmat(diag(D)',n,1));
u=reshape(V*Z*V',n^2,1);
end