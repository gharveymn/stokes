function StokesSFJacobi
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	
	par = Parameters;
	h = par.h;
	toPlot = par.toPlot;
	
	[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,valind,on,xmeshfull,ymeshfull] = MakeGrids;
	
	%note: numel(q) = numel(xmesh) = numel(ymesh)
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
	
	
	nx = numel(xinit);
	ny = numel(yinit);
	
	%make derivative matrices
	bih = biharmonic2(nx,ny,h);
	
	sz = size(bih,1);
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
%	bih = ~(bcw&bce&bcn&bcs).*bih + spdiags(bcw&bce&bcs&bcn,0,sz,sz);
	
	bih = ~onpf.*bih + spdiags(onpf,0,sz,sz);
	
	[L,D,U] = ldu(bih);
	
	Dinv = D^(-1);
	
	disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	vecn = bih\rhs;
	
	q = vecn(1:sz);
	q = filterMat*q;
	
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
	
	u = dy*q;
	v = -dx*q;

	umesh = filterMat'*u;
	Umesh = reshape(umesh,[nx,ny])';

	vmesh = filterMat'*v;
	Vmesh = reshape(vmesh,[nx,ny])';

	qmesh = filterMat'*q;
	Qmesh = reshape(qmesh,[nx,ny])';

	mat = cat(3,Xmesh,Ymesh,Umesh,Vmesh,Qmesh);
	vec = cat(2,xmesh,ymesh,umesh,vmesh,qmesh);

	figs = Plot(mat,vec,par.toPlot,par.filter);
	drawnow;
	
	for i=1:1000
		vecn1 = -Dinv*(L + U)*vecn + Dinv*rhs;
		disp(norm(vecn-vecn1))
		vecn = vecn1;
	
		q = vecn(1:sz);
		q = filterMat*q;

		%make some derivative operator matrices
		%TODO: just make these into a function in the path

		u = dy*q;
		v = -dx*q;

		umesh = filterMat'*u;
		Umesh = reshape(umesh,[nx,ny])';

		vmesh = filterMat'*v;
		Vmesh = reshape(vmesh,[nx,ny])';

		qmesh = filterMat'*q;
		Qmesh = reshape(qmesh,[nx,ny])';

		mat = cat(3,Xmesh,Ymesh,Umesh,Vmesh,Qmesh);
		vec = cat(2,xmesh,ymesh,umesh,vmesh,qmesh);

		Plot(mat,vec,par.toPlot,par.filter,figs);
		drawnow;
	
	end
	
end

