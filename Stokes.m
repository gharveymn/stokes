function Stokes
	%STOKES simulates stokes flow
	
	[x,y,ux,uy] = CreateGrid;
	
	minx = min(x);
	maxx = max(x);
	miny = min(y);
	maxy = max(y);
	absmin = min(minx,miny);
	absmax = max(maxx,maxy);
	
	centerx = (maxx+minx)/2;
	centery = (maxy+miny)/2;
	difa = absmax-absmin;
	
	s = sqrt(ux.^2+uy.^2);
	clrs = s./(max(s));
	cmap = [clrs,zeros(numel(clrs),1),1-clrs];
	
	figure(1)
	%scatter(xcoords,ycoords,[],'.')
	p1 = scatter(x,y,[],cmap,'.');
	axis([centerx-difa/2, centerx+difa/2, centery-difa/2, centery+difa/2])
	
	figure(2)
	p2 = scatter3(x,y,sqrt(ux.^2+uy.^2),[],cmap,'.');
	axis([centerx-difa/2, centerx+difa/2, centery-difa/2, centery+difa/2, -1, 1])	
	
	
end

