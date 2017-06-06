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