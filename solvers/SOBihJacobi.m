function [qmesh,mats] = SOBihJacobi(grids,filtering,rhs,bc,mats)
	
	if(nargin == 7)
		bih = mats{1};
		L = mats{2};
		Dinv = mats{3};
		U = mats{4};
	else
		%make derivative matrices
		bih = biharmonic2(nx,ny,h);
		bih = filterMat*bih*filterMat';
		
		sz = size(bih,1);
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
		
		[L,D,U] = ldu(bih);
		
		Dinv = D^(-1);
		
		mats = {bih,L,Dinv,U};
	end
	
	disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	vecn = bih\rhs;
	
	for i=1:1000
		vecn = -Dinv*(L + U)*vecn + Dinv*rhs;
	end
	qmesh = vecn(1:sz);
end
