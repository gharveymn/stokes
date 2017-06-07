testrun = true;
hs = [0.01;0.02;0.04;0.05;0.1];
Us = {};
xs = zeros(1,numel(hs)-1);

par = Parameters;
par.toPlot = false;
par.ghostpoints = true;
par.filter = true;

par.h = hs(1);
Main
Us{1} = U;
X0 = X;
Y0 = Y;

for i=2:numel(hs)
	par.h = hs(i);
	Main
	Us{i} = U;
	xs(i-1) = getinfnorm(Us{1},Us{i},hs(1),hs(i));
end

plot(hs(2:end),xs);

clear par testrun

function x = getinfnorm(M0,M1,h0,h1)
	%GETINFNORM works for matrices
	%h0 must divide h1
	%ghostpoints are NOT included in the matrices, ie par.filter=true
	n = int32(h1/h0);
	M0 = M0(1:n:end,1:n:end);
	A0 = M0(isfinite(M0));
	A1 = M1(isfinite(M1));
	x = norm(A0-A1,inf);
end