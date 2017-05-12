function laplacian = LaplacianFactory(xsz,ysz,h,neumann)
	% LAPLACIANFACTORY
	%	Produces an nmxnm Laplacian operating on a discrete distribution with spacing h.
	%	This is for form v = [... v_{i,j} v_{i+1,j} v_{i+2,j} ...].
	
	Dxx = sptoeplitz([-2 1],xsz)./h.^2;
	Dyy = sptoeplitz([-2 1],ysz)./h.^2;
	
	if(neumann)
		Dxx(1,1) = -1/h.^2;
		Dxx(1,2) = 1/h.^2;
		Dxx(end,end-1) = -1/h.^2;
		Dxx(end,end) = 1/h.^2;
		
		Dyy(1,1) = -1/h.^2;
		Dyy(1,2) = 1/h.^2;
		Dyy(end,end-1) = -1/h.^2;
		Dyy(end,end) = 1/h.^2;
	end
		
	laplacian = kron(Dyy,speye(xsz))+kron(speye(ysz),Dxx);
	
end

