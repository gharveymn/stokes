function figs = InPost(grids,psimesh,xsz,ysz,filterMat,on,par,figs)
	%INPOST does the post processing of calculation
	
	h = par.h;
	
	%TODO change so that derivatives are cool at boundaries
	%Not strictly necessary for 0 dirichlet boundaries
	
	Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
	Dx(1,:) = 0;
	Dx(end,:) = 0;
	dx = kron(speye(ysz),Dx);
	dx = filterMat*dx*filterMat';
	dx = dx.*~(on);
	
	Dy = sptoeplitz([0 -1],[0 1],ysz)./(2*h);
	Dy(1,:) = 0;
	Dy(end,:) = 0;
	dy = kron(Dy,speye(xsz));
	dy = filterMat*dy*filterMat';
	dy = dy.*~(on);
	
	umesh = dy*psimesh;
	vmesh = -dx*psimesh;
	
	umeshfull = filterMat'*umesh;
	Umesh = reshape(umeshfull,[xsz,ysz])';
	
	vmeshfull = filterMat'*vmesh;
	Vmesh = reshape(vmeshfull,[xsz,ysz])';
	
	psimeshfull = filterMat'*psimesh;
	Psimesh = reshape(psimeshfull,[xsz,ysz])';
	
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

