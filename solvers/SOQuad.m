function [qmesh,mats] = SOQuad(grids,filtering,rhs,bc,mats)
	
	
	if(nargin == 7)
		M = mats{1};
		
	else
		
		%make derivative matrices
		
		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		Dx(1,:) = Dx(2,:);
		Dx(end,:) = Dx(end-1,:);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		
		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		Dy(1,:) = Dy(2,:);
		Dy(end,:) = Dy(end-1,:);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		
		sz = numel(rhs);
		
		Z = sparse(sz,sz);
		I = speye(sz);
		
		M = [-dx -I+dx dy -I+dx I+dx dy
			-dy dx -I+dy -I+dy dx I+dy
			Z Z Z dx -I Z
			Z Z Z dy Z -I
			Z Z Z Z dx dy
			Z dx dy -I Z Z
			dx I Z Z Z Z
			dy Z I Z Z Z];
		
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		M(1:sz,1:sz) = ~(bcinds).*M(1:sz,1:sz) + spdiags(bcinds,0,sz,sz);
		M(1:sz,sz+1:end) = ~(bcinds).*M(1:sz,sz+1:end);
		M(sz+1:2*sz,1:sz) = ~(bcinds).*M(sz+1:2*sz,1:sz) + spdiags(bcinds,0,sz,sz);
		M(sz+1:2*sz,sz+1:end) = ~(bcinds).*M(sz+1:2*sz,sz+1:end);
		
		mats = {M};
		
	end
	
	z = zeros(sz,1);
	rhs = [rhs;rhs;z;z;rhs;z;z;z];
	
	%disp(['lower bound for condition number: ' num2str(condest(M))])
	
	%[L,U] = ilu(M);
	%[vec,flag,relres,iter,resvec] = pcg(M,rhs,1e-8,100,L,U);
	vec = M\rhs;
	qmesh = vec(1:sz);
	
end

