function psimesh = DDASch(grids,filtering,rhs,bcinds,par,solver,h)
	%DDASCH additive schwarz
	
	nx = numel(grids{1});
	ny = numel(grids{2});
	
	psimesh = SODuo(numel(grids{1}),numel(grids{2}),bcinds,rhs,filtering{1},h);
	figs = InPost(grids,psimesh,nx,ny,filtering,par);
	
	%TODO complete the bcinds
	[newgrids,newpsis,newbcinds,newrhss] = Decompose(grids,psimesh,rhs,filtering,par.ddbounds);
	
	p10 = newpsis{1};
	p20 = newpsis{2};
	p30 = newpsis{3};
	
	p11 = 0*newpsis{1};
	p21 = 0*newpsis{2};
	p31 = 0*newpsis{3};
	
	% p21 <- f(p10,p30)
	% p11 <- f(p20)
	% p21 <- f(p30)
	
	nx1 = numel(newgrids{1}{1});
	ny1 = numel(newgrids{1}{2});
	
	nx2 = numel(newgrids{2}{1});
	ny2 = numel(newgrids{2}{2});
	
	nx3 = numel(newgrids{3}{1});
	ny3 = numel(newgrids{3}{2});
	
	%for memory allocation purposes
	rhs1 = newrhss{1};
	rhs2 = newrhss{2};
	rhs3 = newrhss{3};
	
	bcinds0 = newbcinds{1}{1};
	on1 = newbcinds{1}{2};
	onouter1 = newbcinds{1}{3};
	
	bcinds1 = newbcinds{2}{1};
	on2 = newbcinds{2}{2};
	onouter2 = newbcinds{2}{3};
	
	bcinds2 = newbcinds{3}{1};
	%on2 = newbcinds{3}{2};
	%onouter2 = newbcinds{3}{3};
	
	bcinds3 = newbcinds{4}{1};
	on3 = newbcinds{4}{2};
	onouter3 = newbcinds{4}{3};
	
	filterMat1 = speye(nx1*ny1);
	filterMat2 = speye(nx2*ny2);
	filterMat3 = speye(nx3*ny3);
	
	%notation: bcindsx{y} means boundary x (from 0,1,2,3) with respect to the grid y (from 1,2,3)
	
	bcin1 = bcinds1{1} | onouter1;
	bcin2 = bcinds0{2} | bcinds3{2} | onouter2;
	bcin3 = bcinds2{3} | onouter3;
	
	for i = 1:par.dditer
		
		%solve for region 2
		rhs2(bcinds0{2}) = p10(bcinds0{1});
		rhs2(bcinds3{2}) = p30(bcinds3{3});
		p21 = solver(nx2,ny2,bcin2,rhs2,filterMat2,h);
		
		%use solution to get regions 1...
		rhs1(bcinds1{1}) = p21(bcinds1{2});
		p11 = solver(nx1,ny1,bcin1,rhs1,filterMat1,h);
	
		%...and 3
		rhs3(bcinds2{3}) = p21(bcinds2{2});
		p31 = solver(nx3,ny3,bcin3,rhs3,filterMat3,h);
		
		
		%probably put some convergence check here
		p10 = p11;
		p20 = p21;
		p30 = p31;
		
		%now repeat
		psimesh = Recompose(grids,{p11,p21,p31},par.ddbounds);
		figs = InPost(grids,psimesh,nx,ny,filtering,par,figs);
		
		disp(['Iteration ' num2str(i)])
		
	end
	
	psimesh = Recompose(grids,{p11,p21,p31},par.ddbounds);
	
end

