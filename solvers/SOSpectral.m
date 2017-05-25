function psimesh = SOSpectral(xsz,ysz,bcinds,rhs,filterMat,h)

	%make derivative matrices
	
	%change when switch to rectangular
	A = 1/h^2*sptoeplitz([2 -1],xsz);
	%C = sparse([1,xsz],[1,xsz],[2/h^4,2/h^4],xsz,xsz);
	
	AA = sptoeplitz([6/h^4,-4/h^4,1/h^4],xsz);
	[Q,H] = eigs(AA,A,xsz);
	
	G = Q'*reshape(filterMat'*rhs,xsz,xsz)*Q;
	Hii = H*ones(xsz);
	Hjj = Hii';
	
	V = G./(1+Hii+Hjj);
	psi = reshape(P*V*P^-1,xsz*ysz,1);
	psimesh = filterMat*psi;
	
end

