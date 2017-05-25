function psimesh = Iter(xsz,ysz,bcinds,rhs,filterMat,h)
	
	%make derivative matrices
	bih = biharmonic2(xsz,ysz,h);
	bih = filterMat*bih*filterMat';
	
	sz = size(bih,1);
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
	
	Ax = sptoeplitz([2 -1],xsz)./h.^2;
	Ay = sptoeplitz([2 -1],ysz)./h.^2;
	Bx = sptoeplitz([6 -4 1],xsz)./h.^4;
	By = sptoeplitz([6 -4 1],ysz)./h.^4;
	Ix = speye(xsz);
	Iy = speye(ysz);
	
	M = kron(Bx,Iy) + kron(Ix,By) + 2*kron(Ax,Ay);
	M = filterMat*M*filterMat';
	M = ~bcinds.*M + spdiags(bcinds,0,sz,sz);
	
	
	psimesh = pcg(bih,rhs,1e-6,100,M);
	
end

