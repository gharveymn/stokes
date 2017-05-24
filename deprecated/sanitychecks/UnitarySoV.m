function figs = UnitarySoV(figs)
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	
	par = Parameters;
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
	
	%make derivative matrices
	lap = laplacian2(xsz,ysz,h);
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
	psimesh = -vec(1:sz);
	
	if(nargin==1)
		InPost(grids,psimesh,xsz,ysz,filterMat,par,figs);
	else
		figs = InPost(grids,psimesh,xsz,ysz,filterMat,par);
	end
	
end

