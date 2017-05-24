function figs=StokesSFSpectral(figs)
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	par = Parameters;
	h = par.h;
	
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
	H = ymax-ymin;
	bin = 30;
	bout = 0.31654835420315243867796794501877;
	c = -exp(-H^2/(4*bin));
	%func = @(x1,x2) 100*((12.*x1.*exp(-x1.^2./x2))./x2^2 - (8.*x1.^3.*exp(-x1.^2./x2))./x2.^3);
	infunc = @(x1,x2) (2.*x1.*exp(-x1.^2./x2))./x2;
	outfunc = @(x1,x2) (2.*x1.*exp(-x1.^2./x2))./x2;
	
	inflowx = xmin*ones(numel(xmesh),1);
	for i=1:numel(yinit)
		curry = yinit(i);
		ind = (xmesh==inflowx & ymesh==curry);
		in(ind) = infunc(curry,bin);
		bcinds = bcinds | ind;
	end
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	%outflow
	outflowx = xmax*ones(numel(xmesh),1);
	for i=1:numel(yinit)
		curry = yinit(i);
		ind = (xmesh==outflowx & ymesh==curry);
		%out(ind) = 1/30*yinit(i) + 1/12;
		out(ind) = outfunc(curry,bout);
		bcinds = bcinds | ind;
	end
	
	rhs = rhs + out;
	
	for i=1:0
		rhs = rhs + circshift(out,-i);
		bcinds = bcinds | circshift(bcinds,-i);
	end
	
	xsz = numel(xinit);
	ysz = numel(yinit);
	
% 	figure(4)
%  	rmesh = filterMat'*rhs;
%  	Rmesh = reshape(rmesh,[xsz,ysz])';
%  	surf(Xmesh,Ymesh,Rmesh,'edgecolor','none','facecolor','interp');

	%make derivative matrices
	
	%change when switch to rectangular
	A = 1/h^2*sptoeplitz([2 -1],xsz);
	C = sparse([1,xsz],[1,xsz],[2/h^4,2/h^4],xsz,xsz);
	
	
	[P,D] = eigs(C,A,xsz);
	
	G = P'*reshape(filterMat'*rhs,xsz,xsz)*P;
	V = G./(4+repmat(diag(D),1,xsz)+repmat(diag(D)',xsz,1));
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
	
	mat = cat(3,Xmesh,Ymesh,Umesh,Vmesh,Psimesh);
	vec = cat(2,xmeshfull,ymeshfull,umesh,vmesh,psimesh);
	
	if(nargin == 1)
		Plot(mat,vec,par.toPlot,par.filter,figs);
	else
		figs = Plot(mat,vec,par.toPlot,par.filter);
	end
	
end

