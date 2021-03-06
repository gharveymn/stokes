function [qmesh,mats] = SOBihN(grids,filtering,rhs,bc,mats)
	%adds support for neumann boundaries
	
	if(nargin == 7)
		M = mats{1};
	else
		nx = grids{9};
		ny = grids{10};
		h = grids{11};
		filterMat = filtering{1};
		
		%idea: link with dy = 0 on out top/bottom boundaries
		
		%make derivative matrices
		bih = biharmonic2(nx,ny,h,bc{1}{2}{1},bc{1}{2}{2},bc{2}{2}{1},bc{2}{2}{2});
		o = bc{2}{2}{1}&filtering{5}{2};
		bih = spdiag(~o)*bih + circshift(spdiag(o),-1,2);
		%dx = deriv(nx,h,1,2);
		%dy = deriv(ny,h,1,2);
		%qxy = kron(dy,dx);
		%bih = bih + spdiag(bc{2}{2}{1}&~(grids{7}>=5+1.5*grids{11}))*qxy;
		%bih = spdiag(~bc{2}{2}{2})*bih + spdiag(bc{2}{2}{2})*kron(dy,speye(nx));
		M = filterMat*bih*filterMat';
		
		mats = {M};
	end
	
	%disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	qmesh = M\rhs;
	
end

