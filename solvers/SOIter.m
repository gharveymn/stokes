function [psimesh,mats] = SOIter(nx,ny,bcinds,rhs,filterMat,h,mats)
	
	if(nargin == 7)
		
		M = mats{1};
		
	else
		
		%make derivative matrices
		bih = biharmonic2(nx,ny,h);
		bih = filterMat*bih*filterMat';
		
		sz = size(bih,1);
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
		
		Ax = sptoeplitz([2 -1],nx)./h.^2;
		Ay = sptoeplitz([2 -1],ny)./h.^2;
		Bx = sptoeplitz([6 -4 1],nx)./h.^4;
		By = sptoeplitz([6 -4 1],ny)./h.^4;
		Ix = speye(nx);
		Iy = speye(ny);
		
		M = kron(Bx,Iy) + kron(Ix,By) + 2*kron(Ax,Ay);
		M = filterMat*M*filterMat';
		M = ~bcinds.*M + spdiags(bcinds,0,sz,sz);
		
		mats = {M};
		
	end
	
	psimesh = pcg(bih,rhs,1e-6,100,M);
	
end

