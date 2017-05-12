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
	%Lh = ((Lh.*valInd)'.*Lh)';
	
	
	Dx = sptoeplitz([0 -1],[0 1],xsz)./h;
	Dx(1,1) = -1/h;
	Dx(1,2) = 1/h;
	Dy = sptoeplitz([0 -1],[0 1],ysz)./h;
	Dy(1,1) = -1/h;
	Dy(1,2) = 1/h;
	
	ddx = kron(eye(ysz),Dx);
	ddy = kron(eye(xsz),Dy);
	
	%ddx = ((ddx.*valInd)'.*ddx)';
	%ddy = ((ddy.*valInd)'.*ddy)';
	
	fx = zeros(sz,1);
	fy = fx;
	
	for i=1:numel(yinit)
		fx(xmesh==0&ymesh==yinit(i)&valInd(xmesh==0&ymesh==yinit(i))) = 0.5;
	end
	
	zer = sparse(sz,sz);
	
	lhs=[-Lh zer ddx
		zer -Lh ddy
		ddx ddy zer];
	rhs=[fx;fy;zeros(sz,1)];
	
	
	v = lhs\rhs;
	ux = v(1:sz);
	uy = v(sz+1:2*sz);
	p = v(2*sz:end);
	u = sqrt(ux.^2 + uy.^2);
	p = p(valInd);
	
	clrs = p./max(p);
	
	%surf(X,Y,(reshape(u,[xsz,ysz]))','edgecolor','none','facecolor','interp');
	scatter(xmesh(valInd),ymesh(valInd),10,[clrs zeros(numel(p),1) 1-clrs],'.');
	
	minx = min(min(Xmesh));
	maxx = max(max(Xmesh));
	miny = min(min(Ymesh));
	maxy = max(max(Ymesh));
	absmin = min(minx,miny);
	absmax = max(maxx,maxy);
	centerx = (maxx+minx)/2;
	centery = (maxy+miny)/2;
	difa = absmax-absmin;
	axis([centerx-difa/2, centerx+difa/2, centery-difa/2, centery+difa/2])
	
	
end

