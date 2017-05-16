function [xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,valInd,on] = ParseValidIndices
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
		
		xmesh = [];
		ymesh = [];
		
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
			
			xmesh = vertcat(xmesh,currx);
			ymesh = vertcat(ymesh,curry);
			
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
		
		%make grid
		limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
		xinit = (limits(1):h:limits(2))';
		yinit = (limits(3):h:limits(4))';
		xsz = numel(xinit);
		ysz = numel(yinit);
		xmesh = kron(ones(ysz,1),xinit);
		ymesh = kron(yinit,ones(xsz,1));
		
		%Credit to Darren Engwirda for inpoly
		[valInd,on] = inpoly(horzcat(xmesh,ymesh),horzcat(xlimcoords,ylimcoords));
		
		Xmesh = (reshape(xmesh./valInd,[xsz,ysz]))';
		Ymesh = flipud((reshape(ymesh./valInd,[xsz,ysz]))');
		
		xmesh = xmesh(valInd);
		ymesh = ymesh(valInd);
		
	else
		error('Not a valid type of map')
	end
	
end
