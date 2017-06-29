function [qmesh,mats] = SOLap(grids,filtering,rhs,bc,mats)
	
	if(nargin == 7)
		%express lane!
		M = mats{1};
	else
		
		filterMat = filtering{1};
		%make derivative matrices
		lap = laplacian2(grids{9},grids{10},grids{11},[],[],bc{2}{1},bc{2}{2});
		M = filterMat*lap*filterMat';
		
		mats = {M};
		
	end
	
	%disp(['lower bound for condition number: ' num2str(condest(M))])
	
	%[L,U] = ilu(M);
	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	
	qmesh = M\rhs;
	
	
end

