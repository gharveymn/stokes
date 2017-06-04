
if(~exist('par','var'))
	par = Parameters;
end

if(~exist('figs','var'))
	[figs,mat,vec] = run(par);
else
	[figs,mat,vec] = run(par,figs);
end

X = mat(:,:,1);
Y = mat(:,:,2);
U = mat(:,:,3);
V = mat(:,:,4);
Psi = mat(:,:,5);

x = vec(:,1);
y = vec(:,2);
u = vec(:,3);
v = vec(:,4);
psi = vec(:,5);

clear mat vec

function [figs,mat,vec] = run(par,figs)
	
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
	
%  	filterMat = filtering{1};
%   	rmeshfull = filterMat'*rhs;
%   	Rmesh = reshape(rmeshfull,[nx,ny])';
%   	surf(grids{5},grids{6},Rmesh,'edgecolor','none','facecolor','interp');
% 	clr = abs(rmeshfull)./norm(rmeshfull(isfinite(rmeshfull)),inf);
% 	clrs = [clr zeros(numel(clr),1) 1-clr];
%   	scatter3(grids{7},grids{8},rmeshfull,[],clrs,'.');

	if(par.ddrun)
		psimesh = ddsolver(grids,filtering,rhs,bcinds,par,solver,h);
	else
		psimesh = solver(nx,ny,bcinds,rhs,filtering{1},h);
	end
	
	if(exist('figs','var'))
		[figs,mat,vec] = InPost(psimesh,bcinds,grids,filtering,par,figs);
	else
		[figs,mat,vec] = InPost(psimesh,bcinds,grids,filtering,par);
	end
	
end

