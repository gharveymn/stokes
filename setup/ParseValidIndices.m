function [xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,on] = ParseValidIndices
	%PARSEVALIDINDICES parses the map file and gives us our mesh
     %xinit,yinit are unmatched x and y sets -- vector
     %xmesh,ymesh have invalid indices removed -- vector
     %Xmesh,Ymesh have NaN wherever invalid -- matrix
     %filterMat filters out invalids when multiplied, inserts zeros when transposed and multiplied -- matrix
     %on holds the indices which are on the boundary -- in vector form
	
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
     xsz = numel(xinit);
     ysz = numel(yinit);
     xmesh = kron(ones(ysz,1),xinit);
     ymesh = kron(yinit,ones(xsz,1));

     %Credit to Darren Engwirda for inpoly
     [valInd,on] = inpoly(horzcat(xmesh,ymesh),horzcat(xlimcoords,ylimcoords));

     filterMat = spdiags(valInd,0,xsz*ysz,xsz*ysz);
     filterMat = filterMat(valInd,:);
     on = on(valInd);

     Xmesh = (reshape(xmesh./valInd,[xsz,ysz]))';
     Ymesh = flipud((reshape(ymesh./valInd,[ysz,xsz]))');

     xmesh = filterMat*xmesh;
     ymesh = filterMat*ymesh;

	
end
