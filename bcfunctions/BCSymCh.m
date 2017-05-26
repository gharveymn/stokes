function [rhs,bcinds] = BCInOut(xmesh,ymesh,rhs,on)
	%for use with symch map
	bcinds = 0*xmesh;
	
	% add all the the indices which are on the boundary
	bcinds = bcinds | on;
	
	xmax = max(xmesh);
	xmin = min(xmesh);
	
	% inflow
	h = (max(ymesh(xmesh==0))-min(ymesh(xmesh==0)))/2;
	a = 1;
	c = 1/12;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(h^2.*ymesh - ymesh.^3./3) + c;
	in(~(xmesh==inflowx)) = 0;
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	% outflow
	H = max(ymesh(xmesh==xmax))-min(ymesh(xmesh==xmax))/2;
	f = 1/12;
	e = (4*a*h^3/3 + c - f)*3/(4*H^3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(H^2.*ymesh + ymesh.^3./3)+f;
	out(~(xmesh==outflowx)) = 0;
	
	rhs = rhs + out;
	
	for i=1:0
		rhs = rhs + circshift(out,-i);
		bcinds = bcinds | circshift(bcinds,-i);
	end
	
	% set top
	% only for square domain for now
	rhs(ymesh > ymin + H & ~(xmesh==inflowx | xmesh==outflowx)) = max(rhs);
	
end

