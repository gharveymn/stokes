function [grids,filterMat,valind,on] = MakeGrids
	%MakeGrids parses the map file and gives us our mesh
     %xinit,yinit are unmatched x and y sets -- vector
     %xmesh,ymesh have invalid indices removed -- vector
     %Xmesh,Ymesh have NaN wherever invalid -- matrix
     %filterMat filters out invalids when multiplied, inserts zeros when transposed and multiplied -- matrix
     %on holds the indices which are on the boundary -- in vector form, not prefiltered
	
	par = Parameters;
	
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

     %make grid
     limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
     xinit = (limits(1):h:limits(2))';
     yinit = (limits(3):h:limits(4))';
     nx = numel(xinit);
     ny = numel(yinit);
     xmeshfull = kron(ones(ny,1),xinit);
     ymeshfull = kron(yinit,ones(nx,1));

     %Credit to Darren Engwirda for inpoly
     [valind,on] = inpoly(horzcat(xmeshfull,ymeshfull),horzcat(xlimcoords,ylimcoords));

     filterMat = spdiags(valind,0,nx*ny,nx*ny);
     filterMat = filterMat(valind,:);

     Xmesh = (reshape(xmeshfull./valind,[nx,ny]))';
     Ymesh = (reshape(ymeshfull./valind,[ny,nx]))';

     xmesh = filterMat*xmeshfull;
     ymesh = filterMat*ymeshfull;
	
	grids = {xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,xmeshfull,ymeshfull};

	fclose('all');
end
