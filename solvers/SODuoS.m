function [psimesh,mats] = SODuoS(grids,filtering,rhs,bc,mats)
	
	if(nargin == 7)
		%express lane!
		M = mats{1};
	else
		
		
		%make derivative matrices
		lap = laplacian2(nx,ny,h);
		lap = filterMat*lap*filterMat';
		
		sz = size(lap,1);
		
		nw = lap;
		ne = -(speye(sz,sz) + lap);
		sw = sparse(sz,sz);
		se = lap;
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		nw = ~(bcinds).*nw + spdiags(bcinds,0,sz,sz);
		ne = ~(bcinds).*ne;
		
		mats = {M};
		
	end
	
	M = [nw ne
		sw se];
	
	rhs = [rhs;rhs];
	
	%disp(['lower bound for condition number: ' num2str(condest(M))])
	
	%[L,U] = ilu(M);
	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	vec = M\rhs;
	psimesh = vec(1:sz);
	
	
end