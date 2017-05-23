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
	
	
	rmesh = filterMat'*rhs;
	Rmesh = reshape(rmesh,[xsz,ysz])';
	
% 	figure(1)
% 	ax = MakeAxis(Xmesh,Ymesh);
% 	surf(Xmesh,Ymesh,Rmesh,'edgecolor','none','facecolor','interp')
% 	axis(ax)
	
	%make derivative matrices
	lap = laplacian2(xsz,ysz,h);
	
	bcw = 0&on;
	bce = 0&on;
	bcs = 0&on;
	bcn = 0&on;
	bcc = 0&on;
	
	%impose Neumann conditions -- along all boundaries
	%unfortunately I think we need to use loops for this
	
% 	%TODO: think about corners
% 	for i=1:numel(on)
% 		if(on(i))
% 			w = xmeshfull(i)==xmin;
% 			e = xmeshfull(i)==xmax;
% 			s = ymeshfull(i)==ymin;
% 			n = ymeshfull(i)==ymax;
% 			
% 			if(~w)
% 				w = ~valind(i-1);
% 			end
% 			
% 			if(~e)
% 				e = ~valind(i+1);
% 			end
% 			
% 			if(~s)
% 				s = ~valind(i-xsz);
% 			end
% 			
% 			if(~n)
% 				n = ~valind(i+xsz);
% 			end
% 			
% 			%if two of these are satisfied then we have a corner
% 			%cannot be three since Lipschitz, cannot be w&e or n&s for the same reason
% 			%if all false then its on surrounded by all but one side (since it's on the boundary), so corner
% 			
% 			if((w&&n)||(w&&s)||(e&&n)||(e&&s)||~(w||e||s||n))
% 				%corner boundary
% 				bcc(i) = 1;
% 			end
% 			%west boundary
% 			bcw(i) = w;
% 			%east boundary
% 			bce(i) = e;
% 			%south boundary
% 			bcs(i) = s;
% 			%north boundary
% 			bcn(i) = n;
% 		end
% 	end
	
	bcsz = numel(bcw);
	
	%wipe out invalid indices
	lap = filterMat*lap*filterMat';
	bcw = bcw(valind);
	bce = bce(valind);
	bcs = bcs(valind);
	bcn = bcn(valind);
	
	sz = size(lap,1);
	
	nw = lap;
	ne = speye(sz,sz) + lap;
	sw = sparse(sz,sz);
	se = lap;
	
	%impose Dirichlet conditions
	%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
	nw = ~(bcinds|onpf).*nw + spdiags(bcinds|onpf,0,sz,sz);
	ne = ~(bcinds|onpf).*ne;
	se = ~(bcinds|onpf).*se + spdiags(bcinds|onpf,0,sz,sz);
	
	M = [nw ne
		sw se];
	
	rhs = [rhs;rhs];
	
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

