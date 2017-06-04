function [umesh,vmesh,pmesh] = SOPrim(nx,ny,bcinds,rhs,filterMat,h,mats)
	
	fu = rhs{1};
	fv = rhs{2};
	bcu = bcinds{1};
	bcv = bcinds{2};
	
	lap = filterMat*laplacian2(nx,ny,h,4)*filterMat';
	
	%for now. We need to think about the matrix sizes later
	dx = filterMat*deriv(nx*ny,h,1,4)*filterMat';
	dy = filterMat*deriv(nx*ny,h,1,4)*filterMat';
	
	zer = sparse(size(lap,1),size(lap,2));
	
	M11 = -lap;
	M12 = zer;
	M13 = dx;
	M21 = zer;
	M22 = -lap;
	M23 = dy;
	M31 = dx;
	M32 = dy;
	M33 = zer;
	
	M11 = ~bcu.*M11 + spdiag(bcu);
	M13 = ~bcu.*M13;
	M22 = ~bcv.*M22 + spdiag(bcv);
	M23 = ~bcv.*M23;
	
	%maybe M31 and M32 idk
	%M31 = ~bcu.*M31;
	%M32 = ~bcv.*M32 + spdiag(bcv);
	
	
	M = [M11 M12 M13
		M21 M22 M23
		M31 M32 M33];
	
	%L = ichol(M);
	%uvp = pcg(M,[fu;fv;zeros(size(fu))],1e-6,100,L,L');
	uvp = M\[fu;fv;zeros(size(fu))];
	
	umesh = uvp(1:numel(fu));
	vmesh = uvp(numel(fu)+1:2*numel(fu));
	pmesh = uvp(2*numel(fu)+1:end);
	
end

