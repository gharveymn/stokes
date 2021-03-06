function LidCav
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	par = Parameters;
	h = par.h;
	toPlot = par.toPlot;
	
	[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,valind,on,xmeshfull,ymeshfull] = MakeGrids;
	
	%note: numel(q) = numel(xmesh) = numel(ymesh)
	
	xmin = min(xinit);
	xmax = max(xinit);
	ymin = min(yinit);
	ymax = max(yinit);
	
	%make right hand side for Dirichlet BCs
	rhs = 0*xmesh;
	in = 0*xmesh;
	bcinds = rhs;
	
	%driver
	drivery = ymax;
	for i=1:numel(xinit)
		ind = (xmesh == xinit(i) & ymesh==drivery-h);
		in(ind) = 1;
		bcinds = bcinds | ind;
	end
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	
	
	nx = numel(xinit);
	ny = numel(yinit);
	
	%make derivative matrices
	bih = biharmonic2(nx,ny,h);
	
	bcw = 0*on;
	bce = 0*on;
	bcs = 0*on;
	bcn = 0*on;
	bcc = 0*on;
	
	
	%wipe out invalid indices
	bih = filterMat*bih*filterMat';
	
	sz = size(bih,1);
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
	
	on = on(valind);
% 	on(xmesh==0 & ymesh==0.5) = 0;
% 	on(xmesh==0 & ymesh==-0.5) = 0;
% 	on(xmesh==max(xinit) & ymesh==max(yinit)) = 0;
% 	on(xmesh==max(xinit) & ymesh==min(yinit)) = 0;

	bih = ~on.*bih + spdiags(on,0,sz,sz);
	
	q = bih\rhs;
	
	
	%make some derivative operator matrices
	%TODO: just make these into a function in the path
	
	Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
	Dx(1,:) = 0;
	Dx(end,:) = 0;
	dx = kron(speye(ny),Dx);
	dx = filterMat*dx*filterMat';
	dx = dx.*~(sum(dx,2)~=0);
	
	Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
	Dy(1,:) = 0;
	Dy(end,:) = 0;
	dy = kron(Dy,speye(nx));
	dy = filterMat*dy*filterMat';
	dy = dy.*~(sum(dy,2)~=0);
	
	u = -dy*q;
	v = dx*q;
	
	umesh = filterMat'*u;
	Umesh = reshape(umesh,[nx,ny])';
	
	vmesh = filterMat'*v;
	Vmesh = reshape(vmesh,[nx,ny])';
	
	qmesh = filterMat'*q;
	Qmesh = reshape(qmesh,[nx,ny])';
	
	mat = cat(3,Xmesh,Ymesh,Umesh,Vmesh,Qmesh);
	vec = cat(2,xmesh,ymesh,umesh,vmesh,qmesh);
	
	Plot(mat,vec,toPlot);
	
end

