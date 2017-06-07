testrun = true;
hs = [0.01;0.02;0.03;0.04;0.05;0.06;0.1];
Us = zeros(1,numel(hs));
xs = zeros(1,numel(hs)-1);

if(~exist('par','var'))
	par = Parameters;
end

par.h = hs(1);
Main
Us(1) = U;

for i=2:numel(hs)
	par.h = hs(i);
	Main
	xs(i-1) = getinfnorm(Us(1),U,hs(1),hs(i));
end

plot(hs(2:end),xs);

clear par testrun

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