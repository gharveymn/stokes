function [psimesh,mats] = SOLegPseu(grids,filtering,rhs,bc,mats)
	%SOLEGPSEU utilizes a symmetric Legendre pseudospectral method
	
	if(nargin == 7)
		M = mats{1};
		Q = mats{2};
	else
		[xx,rhox,Wx] = lglnodes(nx+1);
		Dx = lpdmat(nx+1,xx);
		Bx = Dx'*Wx*Dx;

		[xy,rhoy,Wy] = lglnodes(ny+1);
		Dy = lpdmat(ny+1,xy);
		By = Dy'*Wy*Dy;

		Px = Bx*Wx^-1*Bx;
		Py = By*Wy^-1*By;
		Px = Px(2:end-1,2:end-1);
		Py = Py(2:end-1,2:end-1);

		Qx = Wx(2:end-1,2:end-1);
		Qy = Wy(2:end-1,2:end-1);
		Q = kron(Qy,Qx);
		Q(logical(diag(bcinds))) = 1;

		Rx = Bx(2:end-1,2:end-1);
		Ry = By(2:end-1,2:end-1);

		M = kron(Py,Qx) + kron(Qy,Px) + 2*kron(Ry,Rx);
		M = ~bcinds.*M + spdiags(bcinds,0,nx*ny,nx*ny);
		
		mats = {M,Q};

	end
	
	G = Q*rhs;
	psimesh = M\G;
	
	
end

