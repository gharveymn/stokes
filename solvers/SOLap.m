function [psimesh,mats] = SOLap(nx,ny,bcinds,rhs,filterMat,h,mats)
	
	if(nargin == 7)
		%express lane!
		M = mats{1};
	else
		
		%make derivative matrices
		lap = -laplacian2(nx,ny,h);
		M = filterMat*lap*filterMat';
		
		sz = size(M,1);
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		M = ~bcinds.*M + spdiags(bcinds,0,sz,sz);
		
		mats = {M};
		
	end
	
	%disp(['lower bound for condition number: ' num2str(condest(M))])
		
	%[L,U] = ilu(M);
	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	
	psimesh = M\rhs;
	
	
end

