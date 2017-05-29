function psimesh = SOSpectral(nx,ny,bcinds,rhs,filterMat,h)

	%make derivative matrices
	%Lui 239,240,241,273,274,281,282
	
	%change when switch to rectangular
	A = 1/h^2*sptoeplitz([2 -1],nx);
	C = sparse([1,nx],[1,nx],[2/h^4,2/h^4],nx,nx);
	A(1,:) = 0;
	A(:,1) = 0;
	A(end,:) = 0;
	A(:,end) = 0;
	A(1,1) = 1;
	A(end,end) = 1;
	
	[Q,H] = eigs(A^2 + C,A,nx);
	
	G = Q'*reshape(filterMat'*rhs,nx,nx)*Q;
	Hii = H*ones(nx);
	Hjj = Hii';
	
	V = G./(1+Hii+Hjj);
	
	
	psi = reshape(Q*V*Q'*A,nx*ny,1);
	psimesh = filterMat*psi;
	
end