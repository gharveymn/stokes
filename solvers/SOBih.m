function [psimesh,mats] = SOBih(grids,filtering,rhs,bc,mats)
	
	if(nargin == 7)
		M = mats{1};
	else
		nx = grids{9};
		ny = grids{10};
		h = grids{11};
		filterMat = filtering{1};
		
		%make derivative matrices
		bih = biharmonic2(nx,ny,h,bc{1}{2}{1},bc{1}{2}{2});
		M = filterMat*bih*filterMat';
		
		mats = {M};
	end
	
	%disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	psimesh = M\rhs;
	
end

