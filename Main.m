if(exist('figs','var'))
	figs = run(figs);
else
	figs = run;
end

function figs = run(figs)
	
	par = Parameters;
	
	rhfunc = par.rhfunc;
	bcfunc = par.bcfunc;
	solver = par.solver;
	h = par.h;
	
	[grids,filtering] = ParseValidIndices(par);
	
	xsz = numel(grids{1});
	ysz = numel(grids{2});
	
	%implement external force function (on rhs)
	rhs = rhfunc(grids{3},grids{4});
	
	%make right hand side for Dirichlet BCs and get indices for those points
	[rhs,bcinds] = bcfunc(grids{3},grids{4},rhs,filtering{3});
	
% 	rmeshfull = filterMat'*rhs;
% 	Rmesh = reshape(rmeshfull,[xsz,ysz])';
% 	surf(grids{5},grids{6},Rmesh,'edgecolor','none','facecolor','interp');
% 	scatter3(grids{7},grids{8},rmeshfull,[],'.');
	
	psimesh = solver(xsz,ysz,bcinds,rhs,filtering{1},h);
	
	if(nargin==1)
		figs = InPost(grids,psimesh,xsz,ysz,filtering,par,figs);
	else
		figs = InPost(grids,psimesh,xsz,ysz,filtering,par);
	end
	
end

