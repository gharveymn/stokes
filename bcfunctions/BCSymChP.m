function [rhs,bc] = BCSymCh(ugrids,filtering,rhs,par)
	
	ugrids = ugrids{3};
	vgrids = ugrids{4};
	
	uxmeshfull = ugrids{7};
	uymeshfull = ugrids{8};
	uxmesh = ugrids{3};
	uymesh = ugrids{4};
	nxp1 = ugrids{9};
	valindinner = filtering{2}{1};
	valindouter = filtering{2}{2};
	
	on = filtering{3}{1};
	onfull = filtering{3}{2};
	dbcfull = filtering{4}{2};
	gpca = filtering{5};
	
	h = par.h;
	
	% add all the the indices which are on the boundary
	% bc{1} - dirichlet
	% bc{2} - neumann
	bc = {{{on,on},{onfull,onfull}},{{[],[]},{[],[]}}};
	
	xmax = max(uxmesh(on));
	xmin = min(uxmesh(on));
	ymax = max(uymesh(on));
	ymin = min(uymesh(on));
	
	% inflow
	inflowmax = max(uymesh(uxmesh==xmin & on));
	inflowmin = min(uymesh(uxmesh==xmin & on));
	
	d = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + d;
	
	a = par.inflowAmp;
	c = 1/12;
	
	inflowx = xmin*ones(numel(uxmesh),1);
	in = a*(d^2.*(uymesh-centerin) - (uymesh-centerin).^3./3) + a*c;
	in(~(uxmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(uymesh(uxmesh==xmax & on));
	outflowmin = min(uymesh(uxmesh==xmax & on));
	
	D = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + D;
	
	f = 1/12;
	e = (4*a*d^3/3 + c - f)*3/(4*D^3);
	
	outflowx = xmax*ones(numel(uxmesh),1);
	out = e*(D^2.*(uymesh-centerout) - (uymesh-centerout).^3./3) + a*f;
	out(~(uxmesh==outflowx) | ~on) = 0;
	
	rhs = rhs + out;
	
	% set top
	rhs(uymesh > centerout & ~(uxmesh==inflowx | uxmesh==outflowx) & on) = out(uxmesh==xmax&uymesh==ymax);
	
	% set bottom
	rhs(uymesh < centerout & ~(uxmesh==inflowx | uxmesh==outflowx) & on) = out(uxmesh==xmax&uymesh==ymin);
	
	bc{1}{1}{1} = bc{1}{1}{1}|gpca{1}(valindouter);
	bc{1}{1}{2} = bc{1}{1}{1};
	bc{1}{2}{1} = logical(filtering{1}'*(1*bc{1}{1}{1}));
	bc{1}{2}{2} = logical(filtering{1}'*(1*bc{1}{1}{2}));
	
	
end

