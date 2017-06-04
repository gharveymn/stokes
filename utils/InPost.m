function [figs,mat,vec] = InPost(psimesh,bcinds,grids,filtering,par,figs)
	%INPOST does the post processing of calculation
	
	h = par.h;
	psimeshfull = filtering{1}'*psimesh;
	
	nx = grids{9};
	ny = grids{10};
	
	bcindsfull = logical(filtering{1}'*(1*bcinds));
	
	%TODO change so that derivatives are cool at boundaries
	if(par.ghostpoints)
		
		if(par.filter)
			[~,~,~,bcindsfull] = closure(grids,filtering,h,'inner',filtering{5}{2},bcindsfull);
			[~,~,~,psimeshfull] = closure(grids,filtering,h,'inner',filtering{5}{2},psimeshfull);
			[~,newgrids,newfiltering,gp] = closure(grids,filtering,h,'inner',filtering{5}{2},filtering{5}{1});
			
		% 	[~,~,~,umeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
		% 	[~,~,~,vmeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
		% 	[~,~,~,psimeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
		% 	[~,newgrids,newfiltering,gp] = closure(newgrids,newfiltering,h,'inner',gp,gp);
		
			[~,~,~,bcindsfull] = closure(newgrids,newfiltering,h,'inner',gp,bcindsfull);
			[~,~,~,psimeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
			[~,grids,filtering,gp] = closure(newgrids,newfiltering,h,'inner',gp,gp);
		end
		
		%TODO figure out how to get back our psi at the right size
		filterMat = filtering{1};
		bc = filtering{4}{1};
		valind = filtering{2}{1};
		on = filtering{3}{1};
		onfull = filtering{3}{2};
		
		bcinds = logical(filterMat*(1*bcindsfull));
		psimesh = filterMat*psimeshfull;
		nx = grids{9};
		ny = grids{10};
		
		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		dx = ~(bc{1}|bc{2}).*dx;
		
		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		dy = ~(bc{3}|bc{4}).*dy;
		
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
	
% 	[~,~,~,umeshfull] = closure(grids,filtering,h,'inner',filtering{5}{2},umeshfull);
% 	[~,~,~,vmeshfull] = closure(grids,filtering,h,'inner',filtering{5}{2},vmeshfull);
% 	[~,~,~,psimeshfull] = closure(grids,filtering,h,'inner',filtering{5}{2},psimeshfull);
% 	[~,newgrids,newfiltering,gp] = closure(grids,filtering,h,'inner',filtering{5}{2},filtering{5}{1});
% 	
% % 	[~,~,~,umeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
% % 	[~,~,~,vmeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
% % 	[~,~,~,psimeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
% % 	[~,newgrids,newfiltering,gp] = closure(newgrids,newfiltering,h,'inner',gp,gp);
% 
% 	[~,~,~,umeshfull] = closure(newgrids,newfiltering,h,'inner',gp,umeshfull);
% 	[~,~,~,vmeshfull] = closure(newgrids,newfiltering,h,'inner',gp,vmeshfull);
% 	[~,~,~,psimeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
% 	[~,grids,filtering,gp] = closure(newgrids,newfiltering,h,'inner',gp,gp);
% 	
% 	filterMat = filtering{1};
% 	valind = filtering{2}{1};
% 	on = filtering{3}{1};
% 	onfull = filtering{3}{2};
% 	psimesh = filterMat*psimeshfull;
% 	nx = grids{9};
% 	ny = grids{10};
% 	
% 	umesh = filterMat*umeshfull;
% 	vmesh = filterMat*vmeshfull;
% 	psimesh = filterMat*psimeshfull;
% 	
% 	Umesh = reshape(umeshfull,[nx,ny])';
% 	Vmesh = reshape(vmeshfull,[nx,ny])';
% 	Psimesh = reshape(psimeshfull,[nx,ny])';
	
	% use matrices rather than cell arrays so they throw a dimension error if we have a bug mismatched
	mat = cat(3,grids{5},grids{6},Umesh,Vmesh,Psimesh);
	vec = cat(2,grids{3},grids{4},umesh,vmesh,psimesh);
	
	if(par.plot)
		if(exist('figs','var'))
			figs = Plot(mat,vec,par,figs);
		else
			figs = Plot(mat,vec,par);
		end
	end
	
end

