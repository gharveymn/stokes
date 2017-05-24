function laplacian = LaplacianFactory(xsz,ysz,h,neumann)
	% LAPLACIANFACTORY
	%	Produces a sparse nmxnm Laplacian operating on a discrete distribution with spacing h.
	%	This is for form v = [... v_{i,j} v_{i+1,j} v_{i+2,j} ...].
     
     %neumann is not neccessary
	
	Dxx = sptoeplitz([-2 1],xsz)./h.^2;
	Dyy = sptoeplitz([-2 1],ysz)./h.^2;
	
	if(nargin == 4 && neumann)
		Dxx(1,:) = Dxx(2,:);
		Dxx(end,:) = Dxx(end-1,:);
		
		Dyy(1,:) = Dyy(2,:);
		Dyy(end,:) = Dyy(end-1,:);
	else
		
	end
	
	laplacian = kron(speye(ysz),Dxx) + kron(Dyy,speye(xsz));
	
end

