A = 1/h^2*sptoeplitz([2 -1],sz);
%C = sparse([1,nx],[1,nx],[2/h^4,2/h^4],nx,nx);

AA = sptoeplitz([6/h^4,-4/h^4,1/h^4],sz);
[Q,H] = eigs(AA,A,sz);
bih = biharmonic2(sz,sz,h);