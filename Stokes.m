function Stokes
	%calculates stokes flow
	
	if (~nargin)
		[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,valInd,on] = ParseValidIndices;
	end
	par = Parameters;
	h = par.h;
	
	sz = numel(valInd);
	xsz = size(Xmesh,2);
	ysz = size(Xmesh,1);
	
	Lh = LaplacianFactory(xsz,ysz,par.h,true);
	%Lh = (Lh'.*valInd)';
	
	Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
	Dx(1,1) = -1/h;
	Dx(1,2) = 1/h;
	Dx(end,end-1) = -1/h;
	Dx(end,end) = 1/h;
	ddx = kron(eye(ysz),Dx);
	
	ddy = sptoeplitz([zeros(1,xsz),-1],[zeros(1,xsz) 1], sz)./(2*h);
	corner = speye(xsz,xsz)./h;
	ddy(1:xsz,1:xsz) = -corner;
	ddy(1:xsz,xsz+1:2*xsz) = corner;
	ddy(end-xsz+1:end,end-2*xsz+1:end-xsz) = -corner;
	ddy(end-xsz+1:end,end-xsz+1:end) = corner;
	
	
	[fx,fy,finds] = GetAppliedForce(xinit,yinit,xmesh,ymesh,valInd,sz);
	
	nLh = -Lh;
	nLh = nLh.*(~on) + spdiags(on,0,size(nLh,1),size(nLh,2));
	%nLh = nLh.*(~finds) + spdiags(finds,0,size(nLh,1),size(nLh,2));
	nLh = nLh(valInd,:)';
	nLh = nLh(valInd,:)';
	
	ddxne = ddx;
	%ddxne = ddx.*(~finds);
	ddxne = ddxne(valInd,:)';
	ddxne = ddxne(valInd,:)';
	
	ddyne = ddy;
	%ddyne = ddy.*(~finds);
	ddyne = ddyne(valInd,:)';
	ddyne = ddyne(valInd,:)';
	
	sz = size(ddxne,1);
	
	zer = sparse(sz,sz);
	
	nw = [nLh zer
		 zer nLh];
	 
	ne = [ddxne;ddyne];
	
	lhs=[nw  ne
		ne' zer];
	
	fx = fx(valInd);
	fy = fy(valInd);
	
	rhs=[fx;fy;zeros(sz,1)];
	
	v = lhs\rhs;
	ux = v(1:sz);
	uy = v(sz+1:2*sz);
	p = v(2*sz+1:end);
	u = sqrt(ux.^2 + uy.^2);
	su = u;
	clrsu = su./max(su);
	
	sp = p;
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
	%surf(Xmesh,Ymesh,(reshape(ux,[xsz,ysz]))','edgecolor','none','facecolor','interp');
	scatter(xmesh(valInd),ymesh(valInd),10,[clrsu zeros(numel(su),1) 1-clrsu],'.');
	axis(ax)
	title('velocity field')
	drawnow
	
	figure(2)
	%surf(Xmesh,Ymesh,(reshape(p,[xsz,ysz]))','edgecolor','none','facecolor','interp');
	scatter(xmesh(valInd),ymesh(valInd),10,[clrsp zeros(numel(sp),1) 1-clrsp],'.');
	axis(ax)
	title('pressure field')
	drawnow
	
	
	
end

