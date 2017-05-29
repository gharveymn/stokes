function Dh = lpdmat(n,x)
	%LPDMAT Summary of this function goes here
	%   Detailed explanation goes here
	
	lnx = legendreP(n,x);
	
	xk = ones(n+1).*x;
	xj = xk';
	lnxk = ones(n+1).*lnx;
	lnxj = lnxk';

	Dh = 1./(xk-xj).*lnxk./lnxj;
	Dh(logical(eye(n+1))) = 0;
	Dh = Dh + diag([(n*(n+1))/4;zeros(n-1,1);-(n*(n+1))/4]);
end

