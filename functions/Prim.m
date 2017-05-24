function psimesh = Prim(xsz,ysz,bcinds,rhs)
	%STOKESSF Calculates Stokes flow using a stream function
	
	%make derivative matrices
	bih = biharmonic2(xsz,ysz,h);
	
	sz = size(bih,1);
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
%	bih = ~(bcw&bce&bcn&bcs).*bih + spdiags(bcw&bce&bcs&bcn,0,sz,sz);
	
	bih = ~onpf.*bih + spdiags(onpf,0,sz,sz);
	
	disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	psimesh = bih\rhs;
	
end

