function [psimesh,M,Q] = SOLegPseu(nx,ny,bcinds,rhs,M,Q)
	%SOLEGPSEU utilizes a symmetric Legendre pseudospectral method
	
	if(nargin >= 7)
		G = Q*rhs;
		psimesh = M\G;
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
		Q = kron(Qx,Qy);
		Q(bcinds,bcinds) = 1;

		Rx = Bx(2:end-1,2:end-1);
		Ry = By(2:end-1,2:end-1);

		G = Q*rhs;

		M = kron(Px,Qy) + kron(Qx,Py) + 2*kron(Rx,Ry);
		M = ~bcinds.*M + spdiags(bcinds,0,sz,sz);

		psimesh = M\G;
	end
	
	
end

