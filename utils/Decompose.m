function [gridsnew,filteringnew,psimeshnew,bcnew,rhsnew,bcb] = Decompose(grids,filtering,psimesh,rhs,par)
	%DECOMPOSE for use with symch
	
	ddbounds = par.ddbounds;
	xmesh = grids{3};
	ymesh = grids{4};
	del = par.h+eps;
	del2 = 2*par.h+eps;
	
	[g1,f1,p1] = getGrid(grids,filtering,psimesh,ddbounds{1},par);
	[g2,f2,p2] = getGrid(grids,filtering,psimesh,ddbounds{2},par);
	[g3,f3,p3] = getGrid(grids,filtering,psimesh,ddbounds{3},par);
	
	gridsnew = {g1,g2,g3};
	filteringnew = {f1,f2,f3};
	psimeshnew = {p1,p2,p3};
	
	on1 = f1{3}{1};
	on2 = f2{3}{1};
	on3 = f3{3}{1};
	
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
	
	%notation: bcx{y} means boundary x located belonging to the region y
	%this is because we need indices for both in each case
	
	
	%bc0 is all with x <= b1 inside R1 on boundary of R2
	bcb0{1} = (xmesh1 >= blx21) & (ymesh1 >= bly21) & (ymesh1 <= bly22) ...
							& ((ymesh1 == bly21) | (ymesh1 == bly22) | (xmesh1 == blx21));
	bcb0{2} = on2 & (xmesh2 <= blx12) & (ymesh2 >= bly11) & (ymesh2 <= bly12);
	
	%belong to R2
	%bc1 is all with x >= b0 inside R2 on boundary of R1
	bcb1{2} = (xmesh2 <= blx12) & (ymesh2 >= bly11) & (ymesh2 <= bly12) ...
						& ((ymesh2 == bly11) | (ymesh2 == bly12) | (xmesh2 == blx12));
	bcb1{1} = on1 & (xmesh1 >= blx21) & (ymesh1 >= bly21) & (ymesh1 <= bly22);
	
	%bc2 is all with x <= b3 inside R2 on boundary of R3
	bcb2{2} = (xmesh2 >= blx31) & (ymesh2 >= bly31) & (ymesh2 <= bly32) ...
							& ((ymesh2 == bly31) | (ymesh2 == bly32) | (xmesh2 == blx31));
	bcb2{3} = on3 & (xmesh3 <= blx22) & (ymesh3 >= bly21) & (ymesh3 <= bly22);
	
	%bc3 is all with x >= b2 inside R3 on boundary of R2
	bcb3{3} = (xmesh3 <= blx22) & (ymesh3 >= bly21) & (ymesh3 <= bly22) ...
						& ((ymesh3 == bly21) | (ymesh3 == bly22) | (xmesh3 == blx22));
	bcb3{2} = on2 & (xmesh2 >= blx31) & (ymesh2 >= bly31) & (ymesh2 <= bly32);
	
	bcsur = filtering{3}{1} | filtering{5}{1}(filtering{2}{2}) | filtering{5}{2}(filtering{2}{2});
	if(par.ghostpoints)
		onouter1 = bcsur((xmesh >= blx11-del2) & (xmesh <= blx12+del2) & (ymesh >= bly11-del2) & (ymesh <= bly12+del2));
		onouter2 = bcsur((xmesh >= blx21-del2) & (xmesh <= blx22+del2) & (ymesh >= bly21-del2) & (ymesh <= bly22+del2));
		onouter3 = bcsur((xmesh >= blx31-del2) & (xmesh <= blx32+del2) & (ymesh >= bly31-del2) & (ymesh <= bly32+del2));
		
		rhs1 = rhs((xmesh >= blx11-del2) & (xmesh <= blx12+del2) & (ymesh >= bly11-del2) & (ymesh <= bly12+del2));
		rhs2 = rhs((xmesh >= blx21-del2) & (xmesh <= blx22+del2) & (ymesh >= bly21-del2) & (ymesh <= bly22+del2));
		rhs3 = rhs((xmesh >= blx31-del2) & (xmesh <= blx32+del2) & (ymesh >= bly31-del2) & (ymesh <= bly32+del2));		
	else
		onouter1 = bcsur((xmesh >= blx11) & (xmesh <= blx12) & (ymesh >= bly11) & (ymesh <= bly12));
		onouter2 = bcsur((xmesh >= blx21) & (xmesh <= blx22) & (ymesh >= bly21) & (ymesh <= bly22));
		onouter3 = bcsur((xmesh >= blx31) & (xmesh <= blx32) & (ymesh >= bly31) & (ymesh <= bly32));
		
		rhs1 = rhs((xmesh >= blx11) & (xmesh <= blx12) & (ymesh >= bly11) &(ymesh <= bly12));
		rhs2 = rhs((xmesh >= blx21) & (xmesh <= blx22) & (ymesh >= bly21) &(ymesh <= bly22));
		rhs3 = rhs((xmesh >= blx31) & (xmesh <= blx32) & (ymesh >= bly31) &(ymesh <= bly32));
	end

	bc1 = bcb1{1} | onouter1;
	bc2 = bcb0{2} | bcb3{2} | onouter2;
	bc3 = bcb2{3} | onouter3;

	%for now
	bc1x = bc1;
	bc1y = bc1;
	bc2x = bc2;
	bc2y = bc2;
	bc3x = bc3;
	bc3y = bc3;
	
	bcnew = {{{bc1x,bc1y},{bc1x,bc1y}},{{bc2x,bc2y},{bc2x,bc2y}},{{bc3x,bc3y},{bc3x,bc3y}}};
	rhsnew = {rhs1,rhs2,rhs3};
	bcb = {bcb0,bcb1,bcb2,bcb3};

end

function [g,f,psimeshnew] = getGrid(grids,filtering,psimesh,bnds,par)
	
	xinit = grids{1};
	yinit = grids{2};
	xmesh = grids{3};
	ymesh = grids{4};
	
	vinds = (xmesh >= bnds{1}(1)-eps) & (xmesh <= bnds{2}(1)+eps) ...
				& (ymesh >= bnds{1}(2)-eps) & (ymesh <= bnds{2}(2)+eps);
	
	%minds = (Xmesh <= bnds{1}(1)) & (Xmesh <= bnds{2}(1)) & (Ymesh <= bnds{1}(2)) & (Ymesh <= bnds{2}(2));
	xinitnew = xinit((xinit >= bnds{1}(1)-eps)&(xinit <= bnds{2}(1)+eps));
	yinitnew = yinit((yinit >= bnds{1}(2)-eps)&(yinit <= bnds{2}(2)+eps));
	
	nxnew = numel(xinitnew);
	nynew = numel(yinitnew);
	
	xmeshnew = xmesh(vinds);
	ymeshnew = ymesh(vinds);
	psimeshnew = psimesh(vinds);
	
	Xmeshnew = reshape(xmeshnew,[nxnew,nynew])';
	Ymeshnew = reshape(ymeshnew,[nxnew,nynew])';
	
	g = {xinitnew,yinitnew,xmeshnew,ymeshnew,Xmeshnew,Ymeshnew,xmeshnew,ymeshnew,nxnew,nynew,grids{11}};

	%filterMat,{valindinner,valindouter},{on,onfull},{dbc,dbcfull},{gp1,gp2}
	
	%these are just rectangles, so we can do this
	filterMatnew = speye(nxnew*nynew);
	valindinnernew = ones(size(xmeshnew));
	valindouternew = valindinnernew;
	onnew = kron([1;zeros(nynew-2,1);1],ones(nxnew,1))|kron(ones(nynew,1),[1;zeros(nxnew-2,1);1]);
	
	f = {filterMatnew,{valindinnernew,valindouternew},{onnew,onnew},{onnew,onnew},{[],[]}};
	%[dbc,dbcfull] = boundarysides(grids,filtering,gp,'outer');

	if(par.ghostpoints)
		[gp1,gnew,fnew] = closure(g,f);
		[~,~,~,psimeshnew] = closure(g,f,[],[],psimeshnew);
		
		[~,~,~,psimeshnew] = closure(gnew,fnew,[],gp1,psimeshnew);
		[gp2,g,f,gp1] = closure(gnew,fnew,[],gp1,gp1);
		f = [f,{gp2}];
		f{5} = {gp1,gp2};
	end
	
end

