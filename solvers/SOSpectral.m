function psimesh = Spectral(xsz,ysz,bcinds,rhs,filterMat,h)

	%make derivative matrices
	
	%change when switch to rectangular
	A = 1/h^2*sptoeplitz([2 -1],xsz);
	C = sparse([1,xsz],[1,xsz],[2/h^4,2/h^4],xsz,xsz);
	
	[P,D] = eigs(A,xsz);
	
	G = P'*reshape(filterMat'*rhs,xsz,xsz)*P;
	Dii = D*ones(xsz);
	Diisq = D^2*ones(xsz);
	Djj = Dii';
	Djjsq = Diisq';
	
		
	
	V = G./(4+repmat(diag(D),1,xsz)+repmat(diag(D)',xsz,1));
	psi = reshape(P*V*P',xsz*ysz,1);
	psimesh = filterMat*psi;
	
end

