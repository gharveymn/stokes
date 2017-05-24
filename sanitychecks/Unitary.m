function Unitary
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	
	par = Parameters;
	h = par.h;
	
	[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,valind,on,xmeshfull,ymeshfull] = ParseValidIndices;
	
	%make right hand side for Dirichlet BCs
	onpf = on(valind);
	rhs = ones(numel(xmesh),1);
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
	psi = -vec(1:sz);
	
	%make some derivative operator matrices
	%TODO: just make these into a function in the path
	
	Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
	Dx(1,:) = 0;
	Dx(end,:) = 0;
	dx = kron(speye(ysz),Dx);
	dx = filterMat*dx*filterMat';
	dx = dx.*~(sum(dx,2)~=0);
	
	Dy = sptoeplitz([0 -1],[0 1],ysz)./(2*h);
	Dy(1,:) = 0;
	Dy(end,:) = 0;
	dy = kron(Dy,speye(xsz));
	dy = filterMat*dy*filterMat';
	dy = dy.*~(sum(dy,2)~=0);
	
	u = dy*psi;
	v = -dx*psi;
	
	umesh = filterMat'*u;
	Umesh = reshape(umesh,[xsz,ysz])';
	
	vmesh = filterMat'*v;
	Vmesh = reshape(vmesh,[xsz,ysz])';
	
	psimesh = filterMat'*psi;
	Psimesh = reshape(psimesh,[xsz,ysz])';
	
	mat = cat(3,Xmesh,Ymesh,Umesh,Vmesh,Psimesh);
	vec = cat(2,xmeshfull,ymeshfull,umesh,vmesh,psimesh);
	
	Plot(mat,vec,par.toPlot,par.filter);
	
end

