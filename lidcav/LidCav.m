function LidCav
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	
	par = Parameters;
	h = par.h;
	toPlot = par.toPlot;
	
	[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,valind,on,xmeshfull,ymeshfull] = ParseValidIndices;
	
	%note: numel(psi) = numel(xmesh) = numel(ymesh)
	
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
	
	inflowx = zeros(numel(xmesh),1);
	for i=1:numel(yinit)
		ind = (xmesh==inflowx & ymesh==yinit(i));
		in(ind) = yinit(i)./4 - yinit(i).^3./3 + 1/12;
		bcinds = bcinds | ind;
	end
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	
	%outflow
	
	outwidth = ymax-ymin;
	outflowx = max(xmesh)*ones(numel(xmesh),1);
	for i=1:numel(yinit)
		ind = (xmesh==outflowx & ymesh==yinit(i));
		out(ind) = 1/(outwidth^3)*(1/4*outwidth^2*yinit(i) - yinit(i)^3/3) + 1/12;
		bcinds = bcinds | ind;
	end
	rhs = rhs + out;
	
	for i=1:0
		rhs = rhs + circshift(out,-i);
		bcinds = bcinds | circshift(bcinds,-i);
	end
	
	
	xsz = numel(xinit);
	ysz = numel(yinit);
	
	
	rmesh = filterMat'*rhs;
	Rmesh = reshape(rmesh,[xsz,ysz])';
	
	figure(1)
	ax = MakeAxis(Xmesh,Ymesh);
	surf(Xmesh,Ymesh,Rmesh,'edgecolor','none','facecolor','interp')
	axis(ax)
	
	%make derivative matrices
	bih = biharmonic2(xsz,ysz,h);
	
	bcw = 0*on;
	bce = 0*on;
	bcs = 0*on;
	bcn = 0*on;
	bcc = 0*on;
	
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
% 			elseif(w)
% 				%west boundary
% 				bcw(i) = 1;
% 			elseif(e)
% 				%east boundary
% 				bce(i) = 1;
% 			elseif(s)
% 				%south boundary
% 				bcs(i) = 1;
% 			elseif(n)
% 				%north boundary
% 				bcn(i) = 1;
% 			else
% 				%debugging
% 				throw(MException('onLoop:boundaryError','entry is on the boundary but not sorted'))
% 			end
% 		end
% 	end
% 	
% 	bcsz = numel(bcw);
% 	o = ones(bcsz,1);
% 	bih = ~bcw.*bih + (o+0.1*rand(bcsz,1)).*bcw.*circshift(bih,-1);
% 	bih = ~bce.*bih + (o+0.1*rand(bcsz,1)).*bce.*circshift(bih,1);
% 	bih = ~bcs.*bih + (o+0.1*rand(bcsz,1)).*bcs.*circshift(bih,-xsz);
% 	bih = ~bcn.*bih + (o+0.1*rand(bcsz,1)).*bcn.*circshift(bih,xsz);
% 	bih = ~bcc.*bih + (o+0.1*rand(bcsz,1)).*spdiags(bcc,0,xsz*ysz,xsz*ysz);
	
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
	
	psi = bih\rhs;
	
	
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

