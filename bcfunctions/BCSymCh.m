function [rhs,bcinds] = BCSymCh(xmesh,ymesh,rhs,on)
	%for use with symch map
	bcinds = 0*xmesh;
	
	% add all the the indices which are on the boundary
	bcinds = bcinds | on;
	
	xmax = max(xmesh);
	xmin = min(xmesh);
	ymax = max(ymesh);
	ymin = min(ymesh);
	
	
	% inflow
	inflowmax = max(ymesh(xmesh==0));
	inflowmin = min(ymesh(xmesh==0));
	h = (inflowmax-inflowmin)/2;
	a = 1;
	c = 1/12;
	
	centerin = inflowmin + h;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(h^2.*(ymesh-centerin) - (ymesh-centerin).^3./3) + c;
	in(~(xmesh==inflowx)) = 0;
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	%outflow
	outflowmax = max(ymesh(xmesh==xmax));
	outflowmin = min(ymesh(xmesh==xmax));
	H = (outflowmax-outflowmin)/2;
	f = 1/12;
	e = (4*a*h^3/3 + c - f)*3/(4*H^3);
	
	centerout = outflowmin + H;
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(H^2.*(ymesh-centerout) - (ymesh-centerout).^3./3) + f;
	out(~(xmesh==outflowx)) = 0;
	
	rhs = rhs + out;
	
	for i=1:0
		rhs = rhs + circshift(out,-i);
		bcinds = bcinds | circshift(bcinds,-i);
	end
	
	% set top
	% only for square domain for now
	rhs(ymesh > ymin + H & ~(xmesh==inflowx | xmesh==outflowx) & on) = max(rhs);
	
end

