if(exist('figs','var'))
	run(figs);
else
	figs = run;
end

function figs = run(figs)
	
	par = Parameters;
	func = par.func;
	h = par.h;
	
	[grids,filterMat,valind,on] = ParseValidIndices;
	
	xinit = grids{1};
	yinit = grids{2};
	
	%make right hand side for Dirichlet BCs
	onpf = on(valind);
	rhs = ones(numel(onpf),1);
	rhs(onpf) = 0;
	
	xsz = numel(xinit);
	ysz = numel(yinit);
	
	bcinds = onpf;
	
	psimesh = func(xsz,ysz,bcinds,rhs,filterMat,h);
	
	if(nargin==1)
		InPost(grids,psimesh,xsz,ysz,filterMat,par,figs);
	else
		figs = InPost(grids,psimesh,xsz,ysz,filterMat,par);
	end
	
end

