function [grids,filtering,par] = ParseValidIndices(par)
	%PARSEVALIDINDICES parses the map file and gives us our mesh
	%xinit,yinit are unmatched x and y sets -- vector
	%xmesh,ymesh have invalid indices removed -- vector
	%Xmesh,Ymesh have NaN wherever invalid -- matrix
	%filterMat filters out invalids when multiplied, inserts zeros when transposed and multiplied -- matrix
	%on holds the indices which are on the boundary -- in vector form, not prefiltered
	
	h = par.h;
	
	file = fopen(par.mapfile, 'r');
	
	formatSpec = '%f';
	data = fscanf(file,formatSpec);
	
	%parse data into coordinates
	xlimcoords = zeros(numel(data)/2,1);
	ylimcoords = xlimcoords;
	
	for i=1:2:numel(data)
		xlimcoords((i+1)/2) = data(i);
		ylimcoords((i+1)/2) = data(i+1);
	end
	
	%make dd bounds (if needed)
	if(par.ddrun)
		par.ddbounds{1}{1}(1) = xlimcoords(1);
		par.ddbounds{1}{2}(1) = xlimcoords(3) + par.ddoverlap;
		par.ddbounds{1}{1}(2) = ylimcoords(1);
		par.ddbounds{1}{2}(2) = ylimcoords(3);
		
		par.ddbounds{2}{1}(1) = xlimcoords(4);
		par.ddbounds{2}{2}(1) = par.h*round(1/par.h*(par.ddmidratio*xlimcoords(5)...
			+(1-par.ddmidratio)*xlimcoords(4))) + par.ddoverlap;
		par.ddbounds{2}{1}(2) = ylimcoords(7);
		par.ddbounds{2}{2}(2) = ylimcoords(4);
		
		par.ddbounds{3}{1}(1) = par.ddbounds{2}{2}(1) - par.ddoverlap;
		par.ddbounds{3}{2}(1) = xlimcoords(5);
		par.ddbounds{3}{1}(2) = ylimcoords(6);
		par.ddbounds{3}{2}(2) = ylimcoords(5);
	end
	
	%make grid
	limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
	xinit = (limits(1):h:limits(2))';
	yinit = (limits(3):h:limits(4))';
	nx = numel(xinit);
	ny = numel(yinit);
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	%Credit to Darren Engwirda for inpoly
	[valind,onfull] = inpoly(horzcat(xmeshfull,ymeshfull),horzcat(xlimcoords,ylimcoords));
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull(valind);
	
	Xmesh = reshape(xmeshfull./valind,[nx,ny])';
	Ymesh = reshape(ymeshfull./valind,[nx,ny])';
	
	xmesh = filterMat*xmeshfull;
	ymesh = filterMat*ymeshfull;
	
	grids = {xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,xmeshfull,ymeshfull,nx,ny};
	
	%filterMat,{valindinner,valindouter},{on,onfull},{bc,bcfull},{gp1,gp2}
	filtering = {filterMat,{valind,valind},{on,onfull},{[],[]},{[],[]}};
	
	if(par.ghostpoints)
		[gp1,grids,filtering] = closure(grids,filtering,h);
		[gp2,grids,filtering,gp1] = closure(grids,filtering,h,'outer',gp1,gp1);
		filtering = [filtering,{gp2}];
	end
	
	%[bc,bcfull] = boundarysides(grids,filtering);
	%filtering{4} = {bc,bcfull};
	filtering{5} = {gp1,gp2};
	
	%we did some arithmetic up there so just
	%make sure ddbounds are actually in the grids
	[~,in] = min(abs(xinit-par.ddbounds{1}{2}(1)));
	par.ddbounds{1}{2}(1) = xinit(in);
	
	[~,in] = min(abs(xinit-par.ddbounds{2}{2}(1)));
	par.ddbounds{2}{2}(1) = xinit(in);
	
	[~,in] = min(abs(xinit-par.ddbounds{3}{1}(1)));
	par.ddbounds{3}{1}(1) = xinit(in);
	
	
	fclose('all');
end
