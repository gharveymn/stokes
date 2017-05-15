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
	Dx(end,:) = Dx(end-1,:);
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
	nLh = nLh.*(~finds) + spdiags(finds,0,size(nLh,1),size(nLh,2));
	
	ddxne = ddx.*(~on) + spdiags(on,0,size(ddx,1),size(ddx,2));
	ddxne = ddxne.*(~finds) + spdiags(finds,0,size(ddxne,1),size(ddxne,2));
	
	ddyne = ddy.*(~on) + spdiags(on,0,size(ddy,1),size(ddy,2));
	ddyne = ddyne.*(~finds) + spdiags(finds,0,size(ddyne,1),size(ddyne,2));
	
	ddxsw = ddx;
	ddysw = ddy;
	
	sz = size(ddxne,1);
	
	zer = sparse(sz,sz);
		
	nw = [nLh zer
		 zer nLh];
	 
	ne = [ddxne;ddyne];
	
	sw = [ddxsw,ddysw];
	
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
	surf(Xmesh,Ymesh,(reshape(ux,[xsz,ysz]))','edgecolor','none','facecolor','interp');
	%scatter(xmesh(valInd),ymesh(valInd),50,[clrsu zeros(numel(su),1) 1-clrsu],'.');
	axis(ax)
	title('velocity field')
	drawnow
	
	figure(2)
	surf(Xmesh,Ymesh,(reshape(p,[xsz,ysz]))','edgecolor','none','facecolor','interp');
	%scatter(xmesh(valInd),ymesh(valInd),50,[clrsp zeros(numel(sp),1) 1-clrsp],'.');
	axis(ax)
	title('pressure field')
	drawnow
	
end

