function StokesSF
	%STOKESSF Calculates Stokes flow using a stream function
	addpath('setup')
	
	[xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,filterMat,on] = ParseValidIndices;
	
	%note: numel(psi) = numel(xmesh) = numel(ymesh)
	
	%make right hand side for Dirichlet BCs
	rhs = 0*xmesh;
	
	%inflow
	bcx = zeros(sz,1);
	bcy = bcx;
	inds = bcy;
	
	for i=1:numel(yinit)
		ind1 = xmesh==zeros(sz,1) & ymesh==yinit(i);
		ind = ind1 & valInd(ind1);
		bcx(ind) = 4*0.3*yinit(i)*(0.5 - yinit(i))/(0.41)^2;
		inds = inds | ind;
	end
	
	
	ax = MakeAxis(Xmesh,Ymesh);
	surf(Xmesh,Ymesh,zeros(size(Xmesh,1),size(Xmesh,2)),'edgecolor','none','facecolor','interp')
	axis(ax)
	
	
end

