if(exist('figs','var'))
	run(figs);
else
	figs = run;
end

function figs = run(figs)
	
	par = Parameters;
	
	rhfunc = par.rhfunc;
	bcfunc = par.bcfunc;
	solver = par.solver;
	h = par.h;
	
	[grids,filterMat,valind,onfull] = ParseValidIndices(par);
	
	xinit = grids{1};
	yinit = grids{2};
	xmesh = grids{3};
	ymesh = grids{4};
	
	xsz = numel(xinit);
	ysz = numel(yinit);
	
	on = onfull(valind);
	
	
	%implement external force function (on rhs)
	rhs = rhfunc(xmesh,ymesh);
	
	%make right hand side for Dirichlet BCs and get indices for those points
	[rhs,bcinds] = bcfunc(xmesh,ymesh,rhs,on);
	
% 	rmeshfull = filterMat'*rhs;
% 	Rmesh = reshape(rmeshfull,[xsz,ysz])';
% 	surf(grids{5},grids{6},Rmesh,'edgecolor','none','facecolor','interp');
% 	scatter3(grids{7},grids{8},rmeshfull,[],'.');
	
	psimesh = solver(xsz,ysz,bcinds,rhs,filterMat,h);
	
	if(nargin==1)
		InPost(grids,psimesh,xsz,ysz,filterMat,par,figs);
	else
		figs = InPost(grids,psimesh,xsz,ysz,filterMat,par);
	end
	
end

