function laplacian = LaplacianFactory(xlen,ylen,h,neumann)
	% LAPLACIANFACTORY
	%	Produces an nmxnm Laplacian operating on a discrete distribution with spacing h.
	
	Ax_h = sptoeplitz([-2 1],xlen)./h.^2;
	Ay_h = sptoeplitz([-2 1],ylen)./h.^2;
	
	if(neumann)
		Ax_h(1,1) = -1/h.^2;
		Ax_h(1,2) = 1/h.^2;
		Ay_h(1,1) = -1/h.^2;
		Ay_h(1,2) = 1/h.^2;
	end
		
	laplacian = kron(eye(ylen),Ax_h)+kron(Ay_h,eye(xlen));
	
end

