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
	
	grids1 = gridsnew{1};
	grids2 = gridsnew{2};
	grids3 = gridsnew{3};
	
	filtering1 = filteringnew{1};
	filtering2 = filteringnew{2};
	filtering3 = filteringnew{3};
	
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
	
	nx1 = grids1{9};
	ny1 = grids1{10};
	
	nx2 = grids2{9};
	ny2 = grids2{10};
	
	nx3 = grids3{9};
	ny3 = grids3{10};
	
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
	rhs2(bc2{1}{1}) = p20(bc2{1}{1});
	rhs2(bcb0{2}) = p10(bcb0{1});
	rhs2(bcb3{2}) = p30(bcb3{3});
	[p21,mats2] = solver(grids2,filtering2,rhs2,bc2);

	%use solution to get regions 1...
	rhs1(bc1{1}{1}) = p10(bc1{1}{1});
	rhs1(bcb1{1}) = p20(bcb1{2});
	[p11,mats1] = solver(grids1,filtering1,rhs1,bc1);

	%...and 3
	rhs3(bc3{1}{1}) = p30(bc3{1}{1});
	rhs3(bcb2{3}) = p20(bcb2{2});
	[p31,mats3] = solver(grids3,filtering3,rhs3,bc3);

	%probably put some convergence check here
	p10 = p11;
	p20 = p21;
	p30 = p31;

	%now repeat
	psimesh = Recompose(grids,{p11,p21,p31},par);
	figs = InPost(psimesh,bc,grids,filtering,par,figs);
	
	for i = 1:par.dditer
		
		if(par.topause > 0)
			pause(par.topause)
		end
		
		%solve for region 2
		rhs2(bc2{1}{1}) = p20(bc2{1}{1});
		rhs2(bcb0{2}) = p10(bcb0{1});
		rhs2(bcb3{2}) = p30(bcb3{3});
		p21 = solver(grids2,filtering2,rhs2,bc2,mats2);
		
		%use solution to get regions 1...
		rhs1(bc1{1}{1}) = p10(bc1{1}{1});
		rhs1(bcb1{1}) = p20(bcb1{2});
		p11 = solver(grids1,filtering1,rhs1,bc1,mats1);
		
		%...and 3
		rhs3(bc3{1}{1}) = p30(bc3{1}{1});
		rhs3(bcb2{3}) = p20(bcb2{2});
		p31 = solver(grids3,filtering3,rhs3,bc3,mats3);
		
		
		%probably put some convergence check here
		p10 = p11;
		p20 = p21;
		p30 = p31;
		
		%now repeat
		psimesh = Recompose(grids,{p11,p21,p31},par);
		figs = InPost(psimesh,bc,grids,filtering,par,figs);
		
		disp(['Iteration ' num2str(i)])
		
		
		
	end
	
	psimesh = Recompose(grids,{p11,p21,p31},par);
	
end

