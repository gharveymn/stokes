function figs=StokesSFDuo(figs)
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	par = Parameters;
	h = par.h;
	
	[grids,filterMat,valind,on] = ParseValidIndices;
	
	%note: numel(psi) = numel(xmesh) = numel(ymesh)
	onpf = on(valind);
	
	xmin = min(xinit);
	xmax = max(xinit);
	ymin = min(yinit);
	ymax = max(yinit);
	
	%make right hand side for Dirichlet BCs
	rhs = 0*xmesh;
	in = 0*xmesh;
	out = 0*xmesh;
	bcinds = rhs;
	
	%inflow
	outwidth = ymax-ymin;
	inflowx = xmin*ones(numel(xmesh),1);
	for i=1:numel(yinit)
		ind = (xmesh==inflowx & ymesh==yinit(i));
		in(ind) = 1/outwidth^3*(outwidth^2/4*yinit(i) - yinit(i)^3/3) + 1/12;
		bcinds = bcinds | ind;
	end
	
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	
	%outflow
	
	outwidth = ymax-ymin;
	outflowx = xmax*ones(numel(xmesh),1);
	for i=1:numel(yinit)
		ind = (xmesh==outflowx & ymesh==yinit(i));
		%out(ind) = 1/30*yinit(i) + 1/12;
		out(ind) = 1/outwidth^3*(outwidth^2/4*yinit(i) - yinit(i)^3/3) + 1/12;
		bcinds = bcinds | ind;
	end
	
	
	rhs = rhs + out;
	
	for i=1:0
		rhs = rhs + circshift(out,-i);
		bcinds = bcinds | circshift(bcinds,-i);
	end
	
	
	xsz = numel(xinit);
	ysz = numel(yinit);
	
	%make derivative matrices
	lap = laplacian2(xsz,ysz,h);
	
	
	sz = size(lap,1);
	
	nw = lap;
	ne = speye(sz,sz) + lap;
	sw = sparse(sz,sz);
	se = lap;
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	nw = ~(bcinds|onpf).*nw + spdiags(bcinds|onpf,0,sz,sz);
	ne = ~(bcinds|onpf).*ne;
	
	M = [nw ne
		sw se];
	
	rhs = [rhs;rhs];
	
	disp(['lower bound for condition number: ' num2str(condest(M))])
	
 	%[L,U] = ilu(M);
 	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	vec = M\rhs;
	psimesh = vec(1:sz);
	
	if(nargin==1)
		InPost(grids,psimesh,xsz,ysz,filterMat,par,figs);
	else
		figs = InPost(grids,psimesh,xsz,ysz,filterMat,par);
	end
	
end

