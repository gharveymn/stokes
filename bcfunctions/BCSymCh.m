function [rhs,bc] = BCSymCh(grids,filtering,rhs,par)
	
	xmeshfull = grids{7};
	ymeshfull = grids{8};
	xmesh = grids{3};
	ymesh = grids{4};
	nx = grids{9};
	valindinner = filtering{2}{1};
	valindouter = filtering{2}{2};
	
	on = filtering{3}{1};
	bcfull = filtering{4}{2};
	gpca = filtering{5};
	
	del = par.h;
	
	%for use with symch map
	bc = 0*xmesh;
	
	% add all the the indices which are on the boundary
	bc = bc | on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	ymax = max(ymesh(on));
	ymin = min(ymesh(on));
	
	
	% inflow
	inflowmax = max(ymesh(xmesh==xmin & on));
	inflowmin = min(ymesh(xmesh==xmin & on));
	
	h = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + h;
	
	
	a = 1;
	c = 1/12;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(h^2.*(ymesh-centerin) - (ymesh-centerin).^3./3) + c;
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(ymesh(xmesh==xmax & on));
	outflowmin = min(ymesh(xmesh==xmax & on));
	
	H = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + H;
	
	f = 1/12;
	e = (4*a*h^3/3 + c - f)*3/(4*H^3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(H^2.*(ymesh-centerout) - (ymesh-centerout).^3./3) + f;
	out(~(xmesh==outflowx) | ~on) = 0;
	
	rhs = rhs + out;
	
	% set top
	rhs(ymesh > centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymax);
	
	% set bottom
	rhs(ymesh < centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymin);
	
	rhs = extendgp(rhs,bcfull,valindouter,gpca,nx);
	bc = bc|gpca{1}(valindouter)|gpca{2}(valindouter);
	
end

