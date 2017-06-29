function [qmesh,mats] = SOBihCG(grids,filtering,rhs,bc,mats)
	
	if(nargin == 7)
		bih = mats{1};
	else
		
		nx = grids{9};
		ny = grids{10};
		h = grids{11};
		filterMat = filtering{1};
		
		%make derivative matrices
		bih = biharmonic2(nx,ny,h);%,bc{2}{1},bc{2}{2});
		bih = filterMat*bih*filterMat';
		
		sz = size(bih,1);
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		%bih = ~bc.*bih + spdiag(bc);
		
		mats = {bih};
	end
	
	%disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	[qmesh,flag,relres,it,resv]=pcg(bih,rhs);
	
end

