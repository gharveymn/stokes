function psimesh = DDMSch(grids,filtering,rhs,bc,par,solver)
	%DDASCH multiplicative schwarz
	
	nx = grids{9};
	ny = grids{10};
	h = grids{11};
	filterMat = filtering{1};
	
	psimesh = SOBih(grids,filtering,rhs,bc);
	figs = InPost(psimesh,bc,grids,filtering,par);
	
	%TODO complete the bc
	[gridsnew,filteringnew,psimeshnew,bcnew,rhsnew,bcb] = Decompose(grids,filtering,psimesh,rhs,par);
	
	p10 = psimeshnew{1};
	p20 = psimeshnew{2};
	p30 = psimeshnew{3};
	
	p11 = 0*psimeshnew{1};
	p21 = 0*psimeshnew{2};
	p31 = 0*psimeshnew{3};
	
	filterMat1 = filteringnew{1};
	filterMat2 = filteringnew{2};
	filterMat3 = filteringnew{3};
	
	
	% p21 <- f(p10,p30)
	% p11 <- f(p20)
	% p21 <- f(p30)
	
	nx1 = gridsnew{1}{9};
	ny1 = gridsnew{1}{10};
	
	nx2 = gridsnew{2}{9};
	ny2 = gridsnew{2}{10};
	
	nx3 = gridsnew{3}{9};
	ny3 = gridsnew{3}{10};
	
	%for memory allocation purposes
	rhs1 = rhsnew{1};
	rhs2 = rhsnew{2};
	rhs3 = rhsnew{3};

	bc1 = bcnew{1};
	bc2 = bcnew{2};
	bc3 = bcnew{3};

	bcb0 = bcb{1};
	bcb1 = bcb{2};
	bcb2 = bcb{3};
	bcb3 = bcb{4};

	%initial calculation
	
	%solve for region 2
	rhs2(bcb0{2}) = p10(bcb0{1});
	rhs2(bcb3{2}) = p30(bcb3{3});
	[p21,mats2] = solver(nx2,ny2,bc2,rhs2,filterMat2,h);

	%use solution to get regions 1...
	rhs1(bcb1{1}) = p21(bcb1{2});
	[p11,mats1] = solver(nx1,ny1,bc1,rhs1,filterMat1,h);

	%...and 3
	rhs3(bcb2{3}) = p21(bcb2{2});
	[p31,mats3] = solver(nx3,ny3,bc3,rhs3,filterMat3,h);

	%probably put some convergence check here
	p10 = p11;
	p20 = p21;
	p30 = p31;

	%now repeat
	psimesh = Recompose(grids,{p11,p21,p31},par.ddbounds);
	figs = InPost(psimesh,bc,grids,filtering,par,figs);
	
	for i = 1:par.dditer
		
		%solve for region 2
		rhs2(bcb0{2}) = p10(bcb0{1});
		rhs2(bcb3{2}) = p30(bcb3{3});
		p21 = solver(nx2,ny2,bcin2,rhs2,filterMat2,h,mats2);
		
		%use solution to get regions 1...
		rhs1(bcb1{1}) = p21(bcb1{2});
		p11 = solver(nx1,ny1,bcin1,rhs1,filterMat1,h,mats1);
	
		%...and 3
		rhs3(bcb2{3}) = p21(bcb2{2});
		p31 = solver(nx3,ny3,bcin3,rhs3,filterMat3,h,mats3);
		
		
		%probably put some convergence check here
		p10 = p11;
		p20 = p21;
		p30 = p31;
		
		%now repeat
		psimesh = Recompose(grids,{p11,p21,p31},par.ddbounds);
		figs = InPost(psimesh,bc,grids,filtering,par,figs);
		
		disp(['Iteration ' num2str(i)])
		
	end
	
	psimesh = Recompose(grids,{p11,p21,p31},par.ddbounds);
	
end

