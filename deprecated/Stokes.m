function Stokes
	%calculates stokes flow
	
	if (~nargin)
		[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,on] = ParseValidIndices;
	end
	
	par = Parameters;
	h = par.h;
	
	sz = size(filterMat,1);
	xsz = size(Xmesh,2);
	ysz = size(Xmesh,1);
	
	Lh = LaplacianFactory(xsz,ysz,par.h,true);
	Lh = (filterMat*(filterMat*Lh)')';
	
	Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
	Dx(1,1) = -1;
	Dx(1,2) = 1;
	Dx(end,end-1) = -1;
	Dx(end,end) = 1;
	ddx = kron(eye(ysz),Dx);
	ddx = (filterMat*(filterMat*ddx)')';
	
	ddy = sptoeplitz([zeros(1,xsz),-1],[zeros(1,xsz) 1], xsz*ysz)./(2*h);
 	corner = speye(xsz,xsz)./h;
 	ddy(1:xsz,1:xsz) = -corner;
 	ddy(1:xsz,xsz+1:2*xsz) = corner;
 	ddy(end-xsz+1:end,end-2*xsz+1:end-xsz) = -corner;
	ddy(end-xsz+1:end,end-xsz+1:end) = corner;
	%ddy(1:xsz,:) = ddy(xsz+1:2*xsz,:);
	%ddy(end-xsz+1:end,:) = ddy(end-2*xsz+1:end-xsz,:);
	ddy = (filterMat*(filterMat*ddy)')';
	
	[fx,fy,finds] = GetAppliedForce(xinit,yinit,xmesh,ymesh,on,sz);
	
	nLh = -Lh;
	nLh = nLh.*(~on) + spdiags(on,0,size(nLh,1),size(nLh,2));
	nLh = nLh.*(~finds) + spdiags(finds,0,size(nLh,1),size(nLh,2));
	
	ddxne = ddx.*(~on).*(~finds) + spdiags(on&finds,0,size(ddx,1),size(ddx,2));
	ddyne = -ddy.*(~on).*(~finds) + spdiags(on&finds,0,size(ddy,1),size(ddy,2));

	ddxsw = ddx;
	ddysw = -ddy;
	
	zer = sparse(size(ddxsw,2),size(ddxne,1));
		
	nw = blkdiag(nLh,-nLh);
	ne = [ddxne;ddyne];
	sw = [ddxsw;ddysw]';
	
	lhs=[nw ne
		sw zer];
	
% 	fx = fx(valInd);
% 	fy = fy(valInd);
	
	rhs=[fx;fy;zeros(sz,1)];
	
	v = lhs\rhs;
	ux = v(1:sz);
	uy = v(sz+1:2*sz);
	p = v(2*sz+1:end);
	u = sqrt(ux.^2 + uy.^2);
	su = u;
	clrsu = su./max(su);
	
	UXmesh = reshape(filterMat'*ux,[xsz,ysz])';
	Pmesh = reshape(filterMat'*p,[xsz,ysz])';
	
	sp = p;
	clrsp = sp./max(sp);
	
	minx = min(min(xmesh));
	maxx = max(max(xmesh));
	miny = min(min(ymesh));
	maxy = max(max(ymesh));
	absmin = min(minx,miny);
	absmax = max(maxx,maxy);
	centerx = (maxx+minx)/2;
	centery = (maxy+miny)/2;
	difa = absmax-absmin;
	ax = [centerx-difa/2, centerx+difa/2, centery-difa/2, centery+difa/2];
	
	figure(1)
	surf(Xmesh,Ymesh,UXmesh,'edgecolor','none','facecolor','interp');
	%scatter(xmesh,ymesh,50,[clrsu zeros(numel(su),1) 1-clrsu],'.');
	axis(ax)
	title('velocity field')
	drawnow
	
	figure(2)
	surf(Xmesh,Ymesh,Pmesh,'edgecolor','none','facecolor','interp');
	%scatter(xmesh,ymesh,50,[clrsp zeros(numel(sp),1) 1-clrsp],'.');
	axis(ax)
	title('pressure field')
	drawnow
	
end

