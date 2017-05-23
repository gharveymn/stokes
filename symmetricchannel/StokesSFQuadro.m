function StokesSFDuo
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
	
% 	figure(1)
% 	ax = MakeAxis(Xmesh,Ymesh);
% 	surf(Xmesh,Ymesh,Rmesh,'edgecolor','none','facecolor','interp')
% 	axis(ax)
	
	%make derivative matrices
	
	Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
	Dx(1,:) = Dx(2,:);
	Dx(end,:) = Dx(end-1,:);
	dx = kron(speye(ysz),Dx);
	
	Dy = sptoeplitz([0 -1],[0 1],ysz)./(2*h);
	Dy(1,:) = Dy(2,:);
	Dy(end,:) = Dy(end-1,:);
	dy = kron(Dy,speye(xsz));
	
	sz = xsz*ysz;
	
	Z = sparse(sz,sz);
	I = speye(sz);
	
	M = [-dx -I+dx dy -I+dx I+dx dy
		-dy dx -I+dy -I+dy dx I+dy
		Z Z Z dx -I Z
		Z Z Z dy Z -I
		Z Z Z Z dx dy
		Z dx dy -I Z Z
		dx I Z Z Z Z
		dy Z I Z Z Z];
	
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	M(1:sz,1:sz) = ~(bcinds|onpf).*M(1:sz,1:sz) + spdiags(bcinds|onpf,0,sz,sz);
	M(1:sz,sz+1:end) = ~(bcinds|onpf).*M(1:sz,sz+1:end);
	M(sz+1:2*sz,1:sz) = ~(bcinds|onpf).*M(sz+1:2*sz,1:sz) + spdiags(bcinds|onpf,0,sz,sz);
	M(sz+1:2*sz,sz+1:end) = ~(bcinds|onpf).*M(sz+1:2*sz,sz+1:end);
	
	z = zeros(sz,1);
	rhs = [rhs;rhs;z;z;rhs;z;z;z];
	
 	%[L,U] = ilu(M);
 	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	vec = M\rhs;
	psi = vec(1:sz);
	
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
	
	Xmesh = Xmesh(3:end-2,3:end-2);
	Ymesh = Ymesh(3:end-2,3:end-2);
	Umesh = Umesh(3:end-2,3:end-2);
	Vmesh = Vmesh(3:end-2,3:end-2);
	Psimesh = Psimesh(3:end-2,3:end-2);
	
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

