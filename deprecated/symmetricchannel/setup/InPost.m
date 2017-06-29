function figs = InPost(grids,qmesh,nx,ny,filterMat,par,figs)
	%INPOST does the post processing of calculation
	
	h = par.h;
	
	Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
	Dx(1,:) = 0;
	Dx(end,:) = 0;
	dx = kron(speye(ny),Dx);
	dx = filterMat*dx*filterMat';
	dx = dx.*~(sum(dx,2)~=0);
	
	Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
	Dy(1,:) = 0;
	Dy(end,:) = 0;
	dy = kron(Dy,speye(nx));
	dy = filterMat*dy*filterMat';
	dy = dy.*~(sum(dy,2)~=0);
	
	umesh = dy*qmesh;
	vmesh = -dx*qmesh;
	
	umeshfull = filterMat'*umesh;
	Umesh = reshape(umeshfull,[nx,ny])';
	
	vmeshfull = filterMat'*vmesh;
	Vmesh = reshape(vmeshfull,[nx,ny])';
	
	qmeshfull = filterMat'*qmesh;
	Qmesh = reshape(qmeshfull,[nx,ny])';
	
	mat = cat(3,grids{5},grids{6},Umesh,Vmesh,Qmesh);
	vec = cat(2,grids{3},grids{4},umesh,vmesh,qmesh);
	
	if(nargin == 1)
		Plot(mat,vec,par,figs);
	else
		figs = Plot(mat,vec,par);
	end
	
end

