function [rhs,bcinds] = BCSymCh(grids,filtering,rhs,par)
	
	xmeshfull = grids{7};
	ymeshfull = grids{8};
	xmesh = grids{3};
	ymesh = grids{4};
	nx = grids{9};
	valindinner = filtering{2}{1};
	valindouter = filtering{2}{2};
	
	on = filtering{3}{1};
	bcfull = filtering{4}{2};
	bcw = bcfull{1};
	bce = bcfull{2};
	bcs = bcfull{3};
	bcn = bcfull{4};
	bcc = bcfull{5};
	
	gp1 = filtering{5}{1};
	gp2 = filtering{5}{2};
	
	del = par.h;
	
	%for use with symch map
	bcinds = 0*xmesh;
	
	% add all the the indices which are on the boundary
	bcinds = bcinds | on;
	
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
	
	%set outer regions
	
	%corners (up to second order)
	rhs = makecorners(rhs,bcc,nx,{gp1,gp2},valindouter,2);	
	
	%west
	bcw1 = circshift(bcw,-1);
	bcw2 = circshift(bcw,-2);
	rhs(bcw1(valindouter)) = rhs(bcw(valindouter));
	rhs(bcw2(valindouter)) = rhs(bcw(valindouter));
	
	%east
	bce1 = circshift(bce,1);
	bce2 = circshift(bce,2);
	rhs(bce1(valindouter)) = rhs(bce(valindouter));
	rhs(bce2(valindouter)) = rhs(bce(valindouter));
	
	%south
	bcs1 = circshift(bcs,-nx);
	bcs2 = circshift(bcs,-2*nx);
	rhs(bcs1(valindouter)) = rhs(bcs(valindouter));
	rhs(bcs2(valindouter)) = rhs(bcs(valindouter));
	
	%north
	bcn1 = circshift(bcn,nx);
	bcn2 = circshift(bcn,2*nx);
	rhs(bcn1(valindouter)) = rhs(bcn(valindouter));
	rhs(bcn2(valindouter)) = rhs(bcn(valindouter));
	
	bcinds = bcinds|gp1(valindouter)|gp2(valindouter);
	
	
end

function rhs = makecorners(rhs,bcc,nx,gpca,valind,m)
	%input i for the order
	
	for k=1:m
		gp = gpca{k};
		for i=-k:k
			for j=-k:k
				bccsw = circshift(circshift(bcc,-i),-j*nx) & gp;
				bccswr = circshift(circshift(bccsw,j*nx),i);

				bccnw = circshift(circshift(bcc,-i),j*nx) & gp;
				bccnwr = circshift(circshift(bccnw,-j*nx),i);

				bccse = circshift(circshift(bcc,i),-j*nx) & gp;
				bccser = circshift(circshift(bccse,j*nx),-i);

				bccne = circshift(circshift(bcc,i),j*nx) & gp;
				bccner = circshift(circshift(bccne,-j*nx),-i);

				rhs(bccsw(valind)) = rhs(bccswr(valind));
				rhs(bccnw(valind)) = rhs(bccnwr(valind));
				rhs(bccse(valind)) = rhs(bccser(valind));
				rhs(bccne(valind)) = rhs(bccner(valind));
			end
		end
	end
end

