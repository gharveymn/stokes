function [psimesh,mats] = SOBihN(grids,filtering,rhs,bc,mats)
	%adds support for neumann boundaries
	
	if(nargin == 7)
		M = mats{1};
	else
		nx = grids{9};
		ny = grids{10};
		h = grids{11};
		filterMat = filtering{1};
		
		%make derivative matrices
		bih = biharmonic2(nx,ny,h,bc{2}{1},bc{2}{2});
		M = filterMat*bih*filterMat';
		
		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		dx = spdiag(bc{3})*dx;
		
		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		dy = spdiag(bc{3})*dy;
		
		M = M + (spdiag(~bc{3})*M + dy);
		
		mats = {M};
	end
	
	%disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	psimesh = M\rhs;
	
end

