function [psimesh,mats] = SODuo(grids,filtering,rhs,bc,mats)
	
	if(nargin == 7)
		%express lane!
		M = mats{1};
	else
		
		filterMat = filtering{1};
		
		%make derivative matrices
		lap = laplacian2(grids{9},grids{10},grids{11},[],[],bc{2}{1},bc{2}{2});
		lap = filterMat*lap*filterMat';
		
		sz = size(lap,1);
		
		nw = lap;
		ne = (speye(sz,sz) + lap);
		sw = sparse(sz,sz);
		se = lap;
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		%nw = ~bcinds.*nw + spdiags(bcinds,0,sz,sz);
		ne = spdiag(~(bc{1}{1}|bc{1}{2}))*ne;
		
		M = [nw ne
			sw se];
		mats = {M};
		
	end
	
	rhs = [rhs;zeros(numel(rhs),1)];
	
	%disp(['lower bound for condition number: ' num2str(condest(M))])
		
	%[L,U] = ilu(M);
	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	
	vec = M\rhs;
	psimesh = vec(1:numel(bc{1}{1}));
	
	
end

