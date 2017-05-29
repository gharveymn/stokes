n = 10;
j = (1:n-2)';
gam = 1/2*sqrt(j.*(j+2)./((j+1/2).*(j+3/2)));
M = diag(gam,1) + diag(gam,-1);
e = eig(M);
x = flipud([-1;e;1]);
lnx = legendreP(n,x);

xk = ones(n+1).*x;
xj = xk';
lnxk = ones(n+1).*lnx;
lnxj = lnxk';

Dh = 1./(xk-xj).*lnxk./lnxj;
Dh(logical(eye(n+1))) = 0;
Dh = Dh + diag([(n*(n+1))/4;zeros(n-1,1);-(n*(n+1))/4]);