function StokesSF
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	par = Parameters;
	h = par.h;
	
	[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,on] = ParseValidIndices;
	
	%note: numel(psi) = numel(xmesh) = numel(ymesh)
	
	%make right hand side for Dirichlet BCs
	rhs = 0*xmesh;
	bcinds = rhs;
	
	%inflow
	
	inflowx = zeros(numel(xmesh),1);
	for i=1:numel(yinit)
		ind = (xmesh==inflowx & ymesh==yinit(i));
		rhs(ind) = yinit(i)./4 - yinit(i).^3./3 + 1/12;
		bcinds = bcinds | ind;
	end
	
	%outflow
	
	outwidth = (max(yinit)-min(yinit));
	outflowx = max(xmesh)*ones(numel(xmesh),1);
	for i=1:numel(yinit)
		ind = (xmesh==outflowx & ymesh==yinit(i));
		rhs(ind) = 1/(outwidth^3)*(1/4*outwidth^2*yinit(i) - yinit(i)^3/3) + 1/12;
		bcinds = bcinds | ind;
	end
	
% 	rhs(xmesh==0 & ymesh==0.5) = 0;
% 	rhs(xmesh==0 & ymesh==-0.5) = 0;
% 	rhs(xmesh==max(xmesh) & ymesh==1.5) = 0;
% 	rhs(xmesh==max(xmesh) & ymesh==-1.5) = 0;
	
	
	xsz = numel(xinit);
	ysz = numel(yinit);
	
	
	rmesh = filterMat'*rhs;
	Rmesh = flipud(reshape(rmesh,[ysz,xsz]));
	
	figure(1)
	ax = MakeAxis(Xmesh,Ymesh);
	surf(Xmesh,Ymesh,Rmesh,'edgecolor','none','facecolor','interp')
	axis(ax)
	
	
	%make derivative matrices
	bih = biharmonic2(xsz,ysz,h);
	
	%wipe out invalid indices
	bih = filterMat*bih*filterMat';
	
	sz = size(bih,1);
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	bih = bih.*(~(bcinds|on)) + spdiags((bcinds|on),0,sz,sz);
	
	psi = bih\rhs;
	
	psimesh = filterMat'*psi;
	Psimesh = flipud(reshape(psimesh,[xsz,ysz])');
	
	figure(1)
	ax = MakeAxis(Xmesh,Ymesh);
	surf(Xmesh,Ymesh,Psimesh,'edgecolor','none','facecolor','interp')
	axis(ax)
	
	
	%make some derivative operator matrices
	%TODO: just make these into a function in the path
		
	Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
	Dx(1,:) = 0;
	Dx(end,:) = 0;
	dx = kron(speye(ysz),Dx);
	dx = filterMat*dx*filterMat';
	
	Dy = sptoeplitz([0 -1],[0 1],ysz)./(2*h);
	Dy(1,:) = 0;
	Dy(end,:) = 0;
	dy = kron(Dy,speye(xsz));
	dy = filterMat*dy*filterMat';
	
	u = dy*psi;
	v = -dx*psi;
	
	umesh = filterMat'*u;
	Umesh = flipud(reshape(umesh,[xsz,ysz])');
	
	vmesh = filterMat'*v;
	Vmesh = flipud(reshape(vmesh,[xsz,ysz])');
	
	figure(1)
	ax = MakeAxis(Xmesh,Ymesh);
	surf(Xmesh,Ymesh,Umesh,'edgecolor','none','facecolor','interp')
	axis(ax)
	
	figure(2)
	ax = MakeAxis(Xmesh,Ymesh);
	surf(Xmesh,Ymesh,Vmesh,'edgecolor','none','facecolor','interp')
	axis(ax)
	
	
end

