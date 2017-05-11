function [x,y,ux,uy] = CreateGrid
	% We'll either use an ascii map file or a connected grid outline (connections counter-clockwise in listed form.
	% Last point should also be the first point).
	
	par = Parameters;
	data = [];
	
	h = par.h;
	
	file = fopen(par.mapfile, 'r');
	
	if(par.maptype(1) == 'a')
		%unfinished; turned out to be way harder than I reckoned
		
		% key:
		%	1		filled region
		%	0		empty region
		%	D		Dirichlet boundary
		%	N		Neumann boundary
		%	(#)>		driven motion governed by #
		
		while ~feof(file)
			data = vertcat(data,fgetl(file));
		end
		data = flipud(data);
		
		x = [];
		y = [];
		
		initflow = zeros(size(data,1),size(data,2));
		
		for i=1:size(data,2)
			currx = [];
			curry = [];
			for j=1:size(data,1)
				curr = data(j,i);
				if(curr == '1' || curr == ')' || curr == '>' || curr == '(' || (i > 1 && data(j,i-1) == '('))
					currx = vertcat(currx,i);
					curry = vertcat(curry,j);
				end
				
				if(i < size(data,2) && data(j,i) == '(')
					initflow(j,i) = str2double(data(j,i+1));
				else
					initflow(j,i) = 0;
				end
			end
			
			x = vertcat(x,currx);
			y = vertcat(y,curry);
			
		end
		
	elseif(par.maptype(1) == 'g')
		
		formatSpec = '%f';
		data = fscanf(file,formatSpec);
		
		%parse data into coordinates	
		xlimcoords = zeros(numel(data)/2,1);
		ylimcoords = xlimcoords;
		
		for i=1:2:numel(data)
			xlimcoords((i+1)/2) = data(i);
			ylimcoords((i+1)/2) = data(i+1);
		end
		
		%turn into vectors
		%xlimvects = diff(xlimcoords);
		%ylimvects = diff(ylimcoords);
		
		%make grid
		limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
		xinit = (limits(1):h:limits(2))';
		yinit = (limits(3):h:limits(4))';
		xmesh = kron(xinit,ones(numel(yinit),1));
		ymesh = kron(ones(numel(xinit),1),yinit);
		
		%rows correspond to each limcoords, columns in those rows are meshpoints
		%xvects = (kron(ones(1,numel(xlimcoords(1:end-1))),xmesh) - kron(ones(numel(xmesh),1),xlimcoords(1:end-1)'))';
		%yvects = (kron(ones(1,numel(ylimcoords(1:end-1))),ymesh) - kron(ones(numel(ymesh),1),ylimcoords(1:end-1)'))';
		
		%if all cross products are positive then point is inside
		%indices = (prod((xlimvects.*yvects - xvects.*ylimvects) > 0))' == 1;
		
		%Credit to Darren Engwirda for inpoly
		indices = inpoly(horzcat(xmesh,ymesh),horzcat(xlimcoords,ylimcoords));
		
		x = xmesh(indices);
		y = ymesh(indices);	
	else
		error('Not a valid type of map')
	end
	
	file = fopen(par.initflowfile, 'r');
	formatSpec = '%f';
	data = fscanf(file,formatSpec);
	
	ux = zeros(numel(x),1);
	uy = ux;
	
	for i=1:4:numel(data)
		ux(x(:)==data(i)&y(:)==data(i+1)) = data(i+2);
		uy(x(:)==data(i)&y(:)==data(i+1)) = data(i+3);
		disp(x(x(:)==data(i)&y(:)==data(i+1)))
		disp(data(i))
		disp(data(i+1))
		disp(data(i+2))
		disp(data(i+3))
		disp(' ')
	end
	
end
