function [psimesh,mats] = SOBih(nx,ny,bcinds,rhs,filterMat,h,mats)
	
	if(nargin == 7)
		bih = mats{1};
	else
		%make derivative matrices
		bih = biharmonic2(nx,ny,h);
		bih = filterMat*bih*filterMat';
		
		sz = size(bih,1);
		
		%impose Dirichlet conditions
		%we do this by just wiping out the row by row multiplication and adding back a diagonal of ones
		bih = ~bcinds.*bih + spdiags(bcinds,0,sz,sz);
		
		mats = {bih};
	end
	
	disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	psimesh = bih\rhs;
	
end

