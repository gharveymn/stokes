function psimesh = SOSpectral(xsz,ysz,bcinds,rhs,filterMat,h)

	%make derivative matrices
	
	%change when switch to rectangular
	A = 1/h^2*sptoeplitz([2 -1],xsz);
	C = sparse([1,xsz],[1,xsz],[2/h^4,2/h^4],xsz,xsz);
	A(1,:) = 0;
	A(:,1) = 0;
	A(end,:) = 0;
	A(:,end) = 0;
	A(1,1) = 1;
	A(end,end) = 1;
	
	[Q,H] = eigs(A^2 + C,A,xsz);
	
	G = Q'*reshape(filterMat'*rhs,xsz,xsz)*Q;
	Hii = H*ones(xsz);
	Hjj = Hii';
	
	V = G./(1+Hii+Hjj);
	psi = reshape(Q*V*Q'*A,xsz*ysz,1);
	psimesh = filterMat*psi;
	
end