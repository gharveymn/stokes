function figs = InPost(grids,psimesh,nx,ny,filtering,par,figs)
	%INPOST does the post processing of calculation
	
	h = par.h;
	
	%TODO change so that derivatives are cool at boundaries
	if(par.ghostpoints)
		
		[~,grids,filtering] = closure(grids{7},grids{8},filtering{4},filtering{2},nx,ny,h,'inner');
		[~,grids,filtering] = closure(grids{7},grids{8},filtering{4},filtering{2},nx,ny,h,'inner');
		
		%TODO figure out how to get back our psi at the right size
		filterMat = filtering{1};
		valind = filtering{2};
		on = filtering{3};
		onfull = filtering{4};
		
		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		dx = ~on.*dx;
		
		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		dy = ~on.*dy;
	else
		filterMat = filtering{1};
		valind = filtering{2};
		on = filtering{3};
		onfull = filtering{4};
		
		%Switch to first order on the boundary
		bc = boundarysides(grids{7},grids{8},onfull,valind,nx);
		bcw = bc{1};
		bce = bc{2};
		bcs = bc{3};
		bcn = bc{4};
		bcc = bc{5};

		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		
		if(par.zeroout)
			dx = ~(bcw|bce|bcc).*dx;
		else
			dx = ~bcw.*dx + 1/h*(-spdiag(bcw) + spdiag(bcw(1:end-1),1));
			dx = ~bce.*dx + 1/h*(-spdiag(bce(2:end),-1) + spdiag(bce));
		end

		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		
		if(par.zeroout)
			dy = ~(bcs|bcn|bcc).*dy;
		else
			dy = ~bcs.*dy + 1/h*(-spdiag(bcs) + spdiag(bcs(1:end-nx),nx));
			dy = ~bcn.*dy + 1/h*(-spdiag(bcn(nx+1:end),-nx) + spdiag(bcn));
		end
	end
	
	umesh = dy*psimesh;
	vmesh = -dx*psimesh;
	
	umeshfull = filterMat'*umesh;
	Umesh = reshape(umeshfull,[nx,ny])';
	
	vmeshfull = filterMat'*vmesh;
	Vmesh = reshape(vmeshfull,[nx,ny])';
	
	psimeshfull = filterMat'*psimesh;
	Psimesh = reshape(psimeshfull,[nx,ny])';
	
	if(par.filter)
		grids{3} = grids{3}(~on);
		grids{4} = grids{4}(~on);
		umesh = umesh(~on);
		vmesh = vmesh(~on);
		psimesh = psimesh(~on);	
	end
	
	mat = cat(3,grids{5},grids{6},Umesh,Vmesh,Psimesh);
	vec = cat(2,grids{3},grids{4},umesh,vmesh,psimesh);
	
	if(exist('figs','var'))
		figs = Plot(mat,vec,par,figs);
	else
		figs = Plot(mat,vec,par);
	end
	
end

