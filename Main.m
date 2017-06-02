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
	ddsolver = par.ddsolver;
	h = par.h;
	
	[grids,filtering,par] = ParseValidIndices(par);
	
	nx = numel(grids{1});
	ny = numel(grids{2});
	
	%implement external force function (on rhs)
	rhs = rhfunc(grids{3},grids{4});
	
	%make right hand side for Dirichlet BCs and get indices for those points
	[rhs,bcinds] = bcfunc(grids,filtering,rhs,par);
	
 	filterMat = filtering{1};
  	rmeshfull = filterMat'*rhs;
  	Rmesh = reshape(rmeshfull,[nx,ny])';
  	surf(grids{5},grids{6},Rmesh,'edgecolor','none','facecolor','interp');
  	scatter3(grids{7},grids{8},rmeshfull,[],'.');

	if(par.ddrun)
		psimesh = ddsolver(grids,filtering,rhs,bcinds,par,solver,h);
	else
		psimesh = solver(nx,ny,bcinds,rhs,filtering{1},h);
	end
	
	if(nargin==1)
		figs = InPost(grids,psimesh,nx,ny,filtering,par,figs);
	else
		figs = InPost(grids,psimesh,nx,ny,filtering,par);
	end
	
end

