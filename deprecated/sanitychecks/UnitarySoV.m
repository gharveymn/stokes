function figs = UnitarySoV(figs)
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	
	par = Parameters;
	h = par.h;
	
	[grids,filterMat,valind,on] = MakeGrids;
	
	xinit = grids{1};
	yinit = grids{2};
	
	%make right hand side for Dirichlet BCs
	onpf = on(valind);
	rhs = ones(numel(onpf),1);
	rhs(onpf) = 0;
	
	nx = numel(xinit);
	ny = numel(yinit);
	
	%make derivative matrices
	lap = laplacian2(nx,ny,h);
	lap = filterMat*lap*filterMat';
	
	sz = size(lap,1);
	
	nw = lap;
	ne = speye(sz,sz) + lap;
	sw = sparse(sz,sz);
	se = lap;
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	nw = ~(onpf).*nw + spdiags(onpf,0,sz,sz);
	ne = ~(onpf).*ne;
	
	M = [nw ne
		sw se];
	
	rhs = [rhs;rhs];
	
	disp(['lower bound for condition number: ' num2str(condest(M))])
	
 	%[L,U] = ilu(M);
 	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	vec = M\rhs;
	qmesh = -vec(1:sz);
	
	if(nargin==1)
		InPost(grids,qmesh,nx,ny,filterMat,par,figs);
	else
		figs = InPost(grids,qmesh,nx,ny,filterMat,par);
	end
	
end

