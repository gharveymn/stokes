function psimesh = SODuoJacobi(xsz,ysz,bcinds,rhs,filterMat,h)
		
	%make derivative matrices
	lap = laplacian2(xsz,ysz,h);
	lap = filterMat*lap*filterMat';
	
	sz = size(lap,1);
	
	nw = lap;
	ne = speye(sz,sz) + lap;
	sw = sparse(sz,sz);
	se = lap;
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	nw = ~(bcinds).*nw + spdiags(bcinds,0,sz,sz);
	ne = ~(bcinds).*ne;
	
	M = [-nw -ne
		sw se];
	
	[L,D,U] = ldu(M);
	
	Dinv = D^(-1);
	
	rhs = [rhs;rhs];
	
	disp(['lower bound for condition number: ' num2str(condest(M))])
	
	vecn = M\rhs;
	
	for i=1:1000
		vecn = -Dinv*(L + U)*vecn + Dinv*rhs;
	end
	
	psimesh = vecn(1:sz);
	
end

