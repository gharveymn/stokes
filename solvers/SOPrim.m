function psimesh = Prim(xsz,ysz,bcinds,rhs,filterMat,h)
	
	%make derivative matrices
	bih = biharmonic2(xsz,ysz,h);
	bih = filterMat*bih*filterMat';
	
	sz = size(bih,1);
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
	
	disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	psimesh = bih\rhs;
	
end

