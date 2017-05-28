function [rhs,bcinds] = BCSymChN(xmesh,ymesh,rhs,on,del,par)
	%for use with symch map with natural smoothing
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
	
	centerin = inflowmin + h;
	
	if(par.ghostpoints)
		h = h - del;
	end
	
	a = .2;
	c = 1/12*a;
	
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
	
	centerout = outflowmin + H;
	outflowx = xmax*ones(numel(xmesh),1);
	
	bcinds = bcinds&(~(xmesh==xmax)|((xmesh==xmax)&(ymesh==outflowmax)|(ymesh==outflowmin)));
	
	rhs(xmesh==xmax&ymesh==outflowmax) = in(xmesh==xmin&ymesh==inflowmax);
	rhs(xmesh==xmax&ymesh==outflowmin) = in(xmesh==xmin&ymesh==inflowmin);
	
	% set top
	rhs(ymesh > centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = in(xmesh==xmin&ymesh==inflowmax);
	
	% set bottom
	rhs(ymesh < centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = in(xmesh==xmin&ymesh==inflowmin);
	
	
	
end

