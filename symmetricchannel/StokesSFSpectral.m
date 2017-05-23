function StokesSF
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	
	par = Parameters;
	h = par.h;
	toPlot = par.toPlot;
	
	[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,valind,on,xmeshfull,ymeshfull] = ParseValidIndices;
	
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
		ind = (xmesh==inflowx+h & ymesh==yinit(i));
		%in(ind) = yinit(i)./4 - yinit(i).^3./3 + 1/12;
		in(ind) = 1/outwidth^3*(outwidth^2/4*yinit(i) - yinit(i)^3/3) + 1/12;
		bcinds = bcinds | ind;
	end
	
	in(onpf) = 0;
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	
	%outflow
	
	outwidth = ymax-ymin;
	outflowx = xmax*ones(numel(xmesh),1);
	for i=1:numel(yinit)
		ind = (xmesh==outflowx-h & ymesh==yinit(i));
		%out(ind) = 1/30*yinit(i) + 1/12;
		out(ind) = 1/outwidth^3*(outwidth^2/4*yinit(i) - yinit(i)^3/3) + 1/12;
		bcinds = bcinds | ind;
	end
	
	out(onpf) = 0;
	
	rhs = rhs + out;
	
	for i=1:0
		rhs = rhs + circshift(out,-i);
		bcinds = bcinds | circshift(bcinds,-i);
	end
	
	
	xsz = numel(xinit);
	ysz = numel(yinit);
	
	%make derivative matrices
	
	%change when switch to rectangular
	A = 1/h^2*sptoeplitz([2 -1],xsz);
	C = sparse([1,xsz],[1,xsz],[(2/h^4),(2/h^4)],xsz,xsz);
	
	[P,D] = eigs(A,xsz);
	
	G = P'*reshape(filterMat'*rhs,xsz,xsz)*P;
	V = G./((4+repmat(diag(D),1,xsz)+repmat(diag(D)',xsz,1)));
	psi = reshape(P*V*P',xsz*ysz,1);
	
	psi = filterMat*psi;
	
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
	
	Xmesh = Xmesh(2:end-1,2:end-1);
	Ymesh = Ymesh(2:end-1,2:end-1);
	Umesh = Umesh(2:end-1,2:end-1);
	Vmesh = Vmesh(2:end-1,2:end-1);
	Psimesh = Psimesh(2:end-1,2:end-1);
	
	if(toPlot == "surf")
		figure(1)
		ax = MakeAxis(Xmesh,Ymesh);
		surf(Xmesh,Ymesh,Umesh,'edgecolor','none','facecolor','interp')
		axis(ax)

		figure(2)
		ax = MakeAxis(Xmesh,Ymesh);
		surf(Xmesh,Ymesh,Vmesh,'edgecolor','none','facecolor','interp')
		axis(ax)

		figure(3)
		ax = MakeAxis(Xmesh,Ymesh);
		surf(Xmesh,Ymesh,Psimesh,'edgecolor','none','facecolor','interp')
		axis(ax)
	else
		figure(1)
		ax = MakeAxis(Xmesh,Ymesh);
		scatter3(xmesh,ymesh,u,10,'.')
		axis(ax)

		figure(2)
		ax = MakeAxis(Xmesh,Ymesh);
		scatter3(xmesh,ymesh,v,10,'.')
		axis(ax)

		figure(3)
		ax = MakeAxis(Xmesh,Ymesh);
		scatter3(xmesh,ymesh,psi,10,'.')
		axis(ax)
	end
	
end

