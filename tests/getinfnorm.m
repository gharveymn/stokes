function x = getinfnorm(M0,M1,h0,h1)
	%GETINFNORM works for matrices
	%h0 must divide h1
	%ghostpoints are NOT included in the matrices, ie par.filter=true
	n = int32(h1/h0);
% 	M0 = M0(3:end-2,3:end-2);
% 	M1 = M1(3:end-2,3:end-2);
	M0 = M0(1:n:end,1:n:end);
	A0 = M0(isfinite(M0));
	A1 = M1(isfinite(M1));
	x = norm(A0-A1,inf);
end