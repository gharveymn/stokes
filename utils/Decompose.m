function [newgrids,newpsis,newbcinds] = Decompose(grids,psimesh,filtering,bounds)
	%DECOMPOSE for use with symch
	
	[g1,p1,on1] = getGrids(grids,psimesh,bounds{1});
	[g2,p2,on2] = getGrids(grids,psimesh,bounds{2});
	[g3,p3,on3] = getGrids(grids,psimesh,bounds{3});
	
	newgrids = {g1,g2,g3};
	newpsis = {p1,p2,p3};
	
	%bounding lines in x direction
	blx11 = bounds{1}{1}(1);		%inlet
	blx12 = bounds{1}{2}(1);		%b1
	bly11 = bounds{1}{1}(2);		%lower y bound
	bly12 = bounds{1}{2}(2);		%upper y bound
	
	
	blx21 = bounds{2}{1}(1);		%b0
	blx22 = bounds{2}{2}(1);		%b3
	bly21 = bounds{2}{1}(2);		%lower y bound
	bly22 = bounds{2}{2}(2);		%upper y bound
	
	blx31 = bounds{3}{1}(1);		%b2
	blx32 = bounds{3}{2}(1);		%outlet
	bly31 = bounds{3}{1}(2);		%lower y bound
	bly32 = bounds{3}{2}(2);		%upper y bound
	
	%overlaps are b0 <= x <= b1 and b2 <= x < b3
	
	xmesh1 = g1{3};
	ymesh1 = g1{4};
	
	xmesh2 = g2{3};
	ymesh2 = g2{4};
	
	xmesh3 = g3{3};
	ymesh3 = g3{4};
	
	%belongs to R1
	%b0 is all with x <= b1 inside R1 on boundary of R2
	bcinds0 = (xmesh1 >= blx21) & (ymesh1 >= bly21) & (ymesh1 <= bly22) ...
							& ((ymesh1 == bly21) | (ymesh1 == bly22) | (xmesh1 == blx21));
	
	%belong to R2
	%b1 is all with x >= b0 inside R2 on boundary of R1
	bcinds1 = (xmesh2 <= blx12) & (ymesh2 >= bly11) & (ymesh2 <= bly12) ...
						& ((ymesh2 == bly11) | (ymesh2 == bly12) | (xmesh2 == blx12));
	
	%b2 is all with x <= b3 inside R2 on boundary of R3
	bcinds2 = (xmesh2 >= blx31) & (ymesh2 >= bly31) & (ymesh2 <= bly32) ...
							& ((ymesh2 == bly31) | (ymesh2 == bly32) | (xmesh2 == blx31));
	
	%b3 is all with x >= b2 inside R3 on boundary of R2
	bcinds3 = (xmesh3 <= blx22) & (ymesh3 >= bly21) & (ymesh3 <= bly22) ...
						& ((ymesh3 == bly21) | (ymesh3 == bly22) | (xmesh3 == blx22));
	
	newbcinds = {{bcinds0,on1},{bcinds1,on2},{bcinds2,on2},{bcinds3,on3}};
	
end

function [g,psimeshnew,onnew] = getGrids(grids,psimesh,bnds)
	
	xinit = grids{1};
	yinit = grids{2};
	xmesh = grids{3};
	ymesh = grids{4};
	
	vinds = (xmesh >= bnds{1}(1)) & (xmesh <= bnds{2}(1)) & (ymesh >= bnds{1}(2)) & (ymesh <= bnds{2}(2));
	
	%minds = (Xmesh <= bnds{1}(1)) & (Xmesh <= bnds{2}(1)) & (Ymesh <= bnds{1}(2)) & (Ymesh <= bnds{2}(2));
	xinitnew = xinit((xinit >= bnds{1}(1))&(xinit <= bnds{2}(1)));
	yinitnew = yinit((yinit >= bnds{1}(2))&(yinit <= bnds{2}(2)));
	
	xsznew = numel(xinitnew);
	ysznew = numel(yinitnew);
	
	xmeshnew = xmesh(vinds);
	ymeshnew = ymesh(vinds);
	psimeshnew = psimesh(vinds);
	
	Xmeshnew = reshape(xmeshnew,[xsznew,ysznew])';
	Ymeshnew = reshape(ymeshnew,[xsznew,ysznew])';
	
	%these are just rectangles, so we can do this
	onnew = kron([1;zeros(ysznew-2,1);1],ones(xsznew,1))|kron(ones(ysznew,1),[1;zeros(xsznew-2,1);1]);
	
	g = {xinitnew,yinitnew,xmeshnew,ymeshnew,Xmeshnew,Ymeshnew,xmeshnew,ymeshnew};
	
end

