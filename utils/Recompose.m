function psimesh = Recompose(grids,newpsis,ddbounds)
	%RECOMPOSE for use with symch
	
	%just to allocate
	psimesh = 0*grids{3};
	
	xmesh = grids{3};
	ymesh = grids{4};
	
	b1 = ddbounds{1};
	b2 = ddbounds{2};
	b3 = ddbounds{3};
	
	%maybe we'll average out overlap, but convergence should take care of that
	psimesh((xmesh >= b1{1}(1)) & (xmesh <= b1{2}(1)) & (ymesh >= b1{1}(2)) & (ymesh <= b1{2}(2))) = newpsis{1};
	psimesh((xmesh >= b2{1}(1)) & (xmesh <= b2{2}(1)) & (ymesh >= b2{1}(2)) & (ymesh <= b2{2}(2))) = newpsis{2};
	psimesh((xmesh >= b3{1}(1)) & (xmesh <= b3{2}(1)) & (ymesh >= b3{1}(2)) & (ymesh <= b3{2}(2))) = newpsis{3};

end

