function [rhs,bcinds] = BCLam(xmesh,ymesh,rhs,on)
	%for use with symch map
	bcinds = 0*xmesh;
	
	%add all the the indices which are on the boundary
	bcinds = bcinds | on;
	
	xmax = max(xmesh);
	xmin = min(xmesh);
	
	%inflow
	inwidth = max(ymesh(xmesh==0))-min(ymesh(xmesh==0));
	inflowx = xmin*ones(numel(xmesh),1);
	in = 1./inwidth.^3.*(inwidth.^2./4.*ymesh - ymesh.^3./3) + 1./12;
	in(~(xmesh==inflowx)) = 0;
	
	rhs = rhs + in;
	
	for i=1:0
		rhs = rhs + circshift(in,i);
		bcinds = bcinds | circshift(bcinds,i);
	end
	
	
	%outflow
	
	outwidth = max(ymesh(xmesh==xmax))-min(ymesh(xmesh==xmax));
	outflowx = xmax*ones(numel(xmesh),1);
	out = 1./outwidth.^3.*(outwidth.^2./4.*ymesh - ymesh.^3./3) + 1./12;
	out(~(xmesh==outflowx)) = 0;
	
	rhs = rhs + out;
	
	for i=1:0
		rhs = rhs + circshift(out,-i);
		bcinds = bcinds | circshift(bcinds,-i);
	end
	
end