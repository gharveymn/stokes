function psimesh = DDMSch(grids,filtering,rhs,bcinds,par,solver)
	%DDMSCH work in progress
	
	h = grids{11};
	ddbounds = par.ddbounds;
	
	%bounding lines
	blx11 = ddbounds{1}{1}(1);		%inlet
	blx12 = ddbounds{1}{2}(1);		%b1
	bly11 = ddbounds{1}{1}(2);		%lower y bound
	bly12 = ddbounds{1}{2}(2);		%upper y bound
	
	blx21 = ddbounds{2}{1}(1);		%b0
	blx22 = ddbounds{2}{2}(1);		%b3
	bly21 = ddbounds{2}{1}(2);		%lower y bound
	bly22 = ddbounds{2}{2}(2);		%upper y bound
	
	blx31 = ddbounds{3}{1}(1);		%b2
	blx32 = ddbounds{3}{2}(1);		%outlet
	bly31 = ddbounds{3}{1}(2);		%lower y bound
	bly32 = ddbounds{3}{2}(2);		%upper y bound
	
	nx = numel(grids{1});
	ny = numel(grids{2});
	sz = numel(bcinds);
	
	lap = -laplacian2(nx,ny,h);
	filterMat = filtering{1};
	lap = filterMat*lap*filterMat';
	lap = ~bcinds.*lap + spdiags(bcinds,0,sz,sz);
	
	xmesh = grids{3};
	ymesh = grids{4};
	
	M = lap.*~(xmesh >= blx11 & xmesh <= blx12 & ymesh >= bly11 & ymesh <= bly12)...
		+ lap.*~(xmesh >= blx21 & xmesh <= blx22 & ymesh >= bly21 & ymesh <= bly22)...
		+ lap.*~(xmesh >= blx31 & xmesh <= blx32 & ymesh >= bly31 & ymesh <= bly32);
	
	psimesh = pcg(lap,rhs,1e-8,100,M);
	
end

