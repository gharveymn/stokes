function [newgrids,newpsis,newbcinds,newrhss] = Decompose(grids,psimesh,rhs,filtering,ddbounds)
	%DECOMPOSE for use with symch
	
	xmesh = grids{3};
	ymesh = grids{4};
	
	[g1,p1,on1] = getGrids(grids,psimesh,ddbounds{1});
	[g2,p2,on2] = getGrids(grids,psimesh,ddbounds{2});
	[g3,p3,on3] = getGrids(grids,psimesh,ddbounds{3});
	
	newgrids = {g1,g2,g3};
	newpsis = {p1,p2,p3};
	
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
	
	%overlaps are b0 <= x <= b1 and b2 <= x < b3
	
	xmesh1 = g1{3};
	ymesh1 = g1{4};
	
	xmesh2 = g2{3};
	ymesh2 = g2{4};
	
	xmesh3 = g3{3};
	ymesh3 = g3{4};
	
	%notation: bcindsx{y} means boundary x located belonging to the region y
	%this is because we need indices for both in each case
	
	
	%b0 is all with x <= b1 inside R1 on boundary of R2
	bcinds0{1} = (xmesh1 >= blx21) & (ymesh1 >= bly21) & (ymesh1 <= bly22) ...
							& ((ymesh1 == bly21) | (ymesh1 == bly22) | (xmesh1 == blx21));
	bcinds0{2} = on2 & (xmesh2 <= blx12) & (ymesh2 >= bly11) & (ymesh2 <= bly12);
	
	%belong to R2
	%b1 is all with x >= b0 inside R2 on boundary of R1
	bcinds1{2} = (xmesh2 <= blx12) & (ymesh2 >= bly11) & (ymesh2 <= bly12) ...
						& ((ymesh2 == bly11) | (ymesh2 == bly12) | (xmesh2 == blx12));
	bcinds1{1} = on1 & (xmesh1 >= blx21) & (ymesh1 >= bly21) & (ymesh1 <= bly22);
	
	%b2 is all with x <= b3 inside R2 on boundary of R3
	bcinds2{2} = (xmesh2 >= blx31) & (ymesh2 >= bly31) & (ymesh2 <= bly32) ...
							& ((ymesh2 == bly31) | (ymesh2 == bly32) | (xmesh2 == blx31));
	bcinds2{3} = on3 & (xmesh3 <= blx22) & (ymesh3 >= bly21) & (ymesh3 <= bly22);
	
	%b3 is all with x >= b2 inside R3 on boundary of R2
	bcinds3{3} = (xmesh3 <= blx22) & (ymesh3 >= bly21) & (ymesh3 <= bly22) ...
						& ((ymesh3 == bly21) | (ymesh3 == bly22) | (xmesh3 == blx22));
	bcinds3{2} = on2 & (xmesh2 >= blx31) & (ymesh2 >= bly31) & (ymesh2 <= bly32);
	
	
	on = filtering{3};
	onouter1 = on((xmesh >= blx11) & (xmesh <= blx12) & (ymesh >= bly11) &(ymesh <= bly12));
	onouter2 = on((xmesh >= blx21) & (xmesh <= blx22) & (ymesh >= bly21) &(ymesh <= bly22));
	onouter3 = on((xmesh >= blx31) & (xmesh <= blx32) & (ymesh >= bly31) &(ymesh <= bly32));
	
	rhs1 = rhs((xmesh >= blx11) & (xmesh <= blx12) & (ymesh >= bly11) &(ymesh <= bly12));
	rhs2 = rhs((xmesh >= blx21) & (xmesh <= blx22) & (ymesh >= bly21) &(ymesh <= bly22));
	rhs3 = rhs((xmesh >= blx31) & (xmesh <= blx32) & (ymesh >= bly31) &(ymesh <= bly32));

	newbcinds = {{bcinds0,on1,onouter1},{bcinds1,on2,onouter2},{bcinds2,on2,onouter2},{bcinds3,on3,onouter3}};
	newrhss = {rhs1,rhs2,rhs3};
	
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
	
	nxnew = numel(xinitnew);
	nynew = numel(yinitnew);
	
	xmeshnew = xmesh(vinds);
	ymeshnew = ymesh(vinds);
	psimeshnew = psimesh(vinds);
	
	Xmeshnew = reshape(xmeshnew,[nxnew,nynew])';
	Ymeshnew = reshape(ymeshnew,[nxnew,nynew])';
	
	%these are just rectangles, so we can do this
	onnew = kron([1;zeros(nynew-2,1);1],ones(nxnew,1))|kron(ones(nynew,1),[1;zeros(nxnew-2,1);1]);
	
	g = {xinitnew,yinitnew,xmeshnew,ymeshnew,Xmeshnew,Ymeshnew,xmeshnew,ymeshnew};
	
end

