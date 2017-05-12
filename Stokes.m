function Stokes
	%calculates stokes flow
	
	if (~nargin)
		[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,valInd] = ParseValidIndices;
	end
	par = Parameters;
	h = par.h;
	
	sz = numel(valInd);
	xsz = size(Xmesh,2);
	ysz = size(Xmesh,1);
	
	
	Lh = LaplacianFactory(xsz,ysz,par.h,true);
	Lh = (Lh'.*valInd)';
	
	
	Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
	Dx(1,1) = -1/h;
	Dx(1,2) = 1/h;
	Dx(end,end-1) = -1/h;
	Dx(end,end) = 1/h;
	Dy = sptoeplitz([0 -1],[0 1],ysz)./(2*h);
	Dy(1,1) = -1/h;
	Dy(1,2) = 1/h;
	Dy(end,end-1) = -1/h;
	Dy(end,end) = 1/h;
	
	ddx = kron(eye(ysz),Dx);
	ddy = kron(eye(xsz),Dy);
	
	ddy = sptoeplitz([zeros(xsz,1);-1],[zeros(1,xsz) 1], sz)./(2*h);
	corner = speye(xsz,xsz)./h;
	ddy(1:xsz,1:xsz) = -corner;
	ddy(1:xsz,xsz+1:2*xsz) = corner;
	ddy(end-xsz+1:end,end-2*xsz+1:end-xsz) = -corner;
	ddy(end-xsz+1:end,end-xsz+1:end) = corner;
	
	
	ddx = (ddx'.*valInd)';
	ddy = (ddy'.*valInd)';
	
	[fx,fy] = GetAppliedForce(xinit,yinit,xmesh,ymesh,valInd,sz);
	
	zer = sparse(sz,sz);
	
	lhs=[-Lh zer ddx
		zer -Lh ddy
		ddx ddy zer];
	
	
	rhs=[fx;fy;zeros(sz,1)];
	
	v = lhs\rhs;
	ux = v(1:sz);
	uy = v(sz+1:2*sz);
	p = v(2*sz+1:end);
	u = sqrt(ux.^2 + uy.^2);
	su = u(valInd);
	clrsu = su./max(su);
	
	sp = p(valInd);
	clrsp = sp./max(sp);
	
	
	
	minx = min(min(Xmesh));
	maxx = max(max(Xmesh));
	miny = min(min(Ymesh));
	maxy = max(max(Ymesh));
	absmin = min(minx,miny);
	absmax = max(maxx,maxy);
	centerx = (maxx+minx)/2;
	centery = (maxy+miny)/2;
	difa = absmax-absmin;
	ax = [centerx-difa/2, centerx+difa/2, centery-difa/2, centery+difa/2];
	
	figure(1)	
	surf(Xmesh,Ymesh,(reshape(u,[xsz,ysz]))','edgecolor','none','facecolor','interp');
	%scatter(xmesh(valInd),ymesh(valInd),10,[clrsu zeros(numel(su),1) 1-clrsu],'.');
	axis(ax)
	title('velocity field')
	drawnow
	
	figure(2)
	surf(Xmesh,Ymesh,(reshape(p,[xsz,ysz]))','edgecolor','none','facecolor','interp');
	%scatter(xmesh(valInd),ymesh(valInd),10,[clrsp zeros(numel(sp),1) 1-clrsp],'.');
	axis(ax)
	title('pressure field')
	drawnow
	
	
	
end

