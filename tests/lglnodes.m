n = 11;
j = (1:n-2)';
gam = 1/2*sqrt(j.*(j+2)./((j+1/2).*(j+3/2)));
M = diag(gam,1) + diag(gam,-1);
e = eig(M);


%FEx
N = 11;

N1=N+1;

x=cos(pi*(0:N)/N)';
P=zeros(N1,N1);

xold=2;

while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
             
end

w=2./(N*N1*P(:,N1).^2);