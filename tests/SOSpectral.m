A = 1/h^2*sptoeplitz([2 -1],sz);
%C = sparse([1,xsz],[1,xsz],[2/h^4,2/h^4],xsz,xsz);

AA = sptoeplitz([6/h^4,-4/h^4,1/h^4],sz);
[Q,H] = eigs(AA,A,sz);
	
