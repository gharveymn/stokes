function [newgrids,newpsis] = Decompose(grids,psimesh,bounds)
	%DECOMPOSE for use with symch
	
	
	[g1,p1] = getGrids(grids,psimesh,bounds{1});
	[g2,p2] = getGrids(grids,psimesh,bounds{2});
	[g3,p3] = getGrids(grids,psimesh,bounds{3});
	
	newgrids = {g1,g2,g3};
	newpsis = {p1,p2,p3};
	
end

function [g,psimeshnew] = getGrids(grids,psimesh,bnds)
	
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
	
	g = {xinitnew,yinitnew,xmeshnew,ymeshnew,Xmeshnew,Ymeshnew,xmeshnew,ymeshnew};
	
end

