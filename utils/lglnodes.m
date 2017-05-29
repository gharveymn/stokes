function [x,rho,W] = lglnodes(n)
	%LGLNODES calculates the n+1 Legendre Gauss-Lobatto nodes for n, also includes weights
	
	%we only need to calculate half of these
	j = (1:n-2)';
	gam = 1/2*sqrt(j.*(j+2)./((j+1/2).*(j+3/2)));
	M = spdiag(gam,1) + spdiag(gam,-1);
	opts.issym = true;
	if(mod(n,2) == 0)
		i = round((n-2)/2);
		
		e = eigs(M,i,'la',opts);
		x = [1;e;0;-flipud(e);-1];
		
		r = 2./(n*(n+1).*legendreP(n,x(1:i+2)).^2);
		rho = [r;flipud(r(1:end-1))];
	else
		i = round((n-1)/2);
		e = eigs(M,i,'la',opts);
		x = [1;e;-flipud(e);-1];
		
		r = 2./(n*(n+1).*legendreP(n,x(1:i+1)).^2);
		rho = [r;flipud(r)];
	end
	
	W = diag(rho);
	
end

