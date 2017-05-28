function figs = InPost(grids,psimesh,xsz,ysz,filtering,par,figs)
	%INPOST does the post processing of calculation
	
	filterMat = filtering{1};
	valind = filtering{2};
	on = filtering{3};
	onfull = filtering{4};
	
	h = par.h;
	
	%TODO change so that derivatives are cool at boundaries
	if(par.ghostpoints)
		Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
		dx = kron(speye(ysz),Dx);
		dx = filterMat*dx*filterMat';
		dx = ~on.*dx;
		
		Dy = sptoeplitz([0 -1],[0 1],ysz)./(2*h);
		dy = kron(Dy,speye(xsz));
		dy = filterMat*dy*filterMat';
		dy = ~on.*dy;
	else
		%Switch to first order on the boundary
		[bcw,bce,bcs,bcn,bcc] = getWhereBoundaries(grids{7},grids{8},onfull,valind,xsz);

		Dx = sptoeplitz([0 -1],[0 1],xsz)./(2*h);
		dx = kron(speye(ysz),Dx);
		dx = filterMat*dx*filterMat';
		dx = ~bcw.*dx + 1/h*(-spdiag(bcw) + spdiag(bcw(1:end-1),1));
		dx = ~bce.*dx + 1/h*(-spdiag(bce(2:end),-1) + spdiag(bce));

		Dy = sptoeplitz([0 -1],[0 1],ysz)./(2*h);
		dy = kron(Dy,speye(xsz));
		dy = filterMat*dy*filterMat';
		dy = ~bcs.*dy + 1/h*(-spdiag(bcs) + spdiag(bcs(1:end-xsz),xsz));
		dy = ~bcn.*dy + 1/h*(-spdiag(bcn(xsz+1:end),-xsz) + spdiag(bcn));
	end
	
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
		Plot(mat,vec,par,figs);
	else
		figs = Plot(mat,vec,par);
	end
	
end

function [bcw,bce,bcs,bcn,bcc] = getWhereBoundaries(xmeshfull,ymeshfull,onfull,valind,xsz)
	
	xmin = min(xmeshfull);
	xmax = max(xmeshfull);
	ymin = min(ymeshfull);
	ymax = max(ymeshfull);
	
	bcw = onfull;
	bce = onfull;
	bcs = onfull;
	bcn = onfull;
	
	xminb = (xmeshfull==xmin);
	xmaxb = (xmeshfull==xmax);
	yminb = (ymeshfull==ymin);
	ymaxb = (ymeshfull==ymax);
	
	bcw = bcw&(xminb);
	bce = bce&(xmaxb);
	bcs = bcs&(yminb);
	bcn = bcn&(ymaxb);
	
	r = circshift(valind&~xmaxb,1);
	l = circshift(valind&~xminb,-1);
	u = circshift(valind&~ymaxb,xsz);
	d = circshift(valind&~yminb,-xsz);
	
	bcw = bcw|onfull&~r;
	bce = bce|onfull&~l;
	bcs = bcs|onfull&~u;
	bcn = bcn|onfull&~d;
	
	%corners--is in two of the previous or is surrounded
	bcc = (bcw&bce)|(bcw&bcs)|(bcw&bcn)|(bce&bcs)|(bce&bcn)|(bcs&bcn);
	
	%inner corner boundary condition
	bcci = onfull&(r&l&u&d);
	
	bcw = bcw|bcci;
	bce = bce|bcci;
	bcs = bcs|bcci;
	bcn = bcn|bcci;
	bcc = bcc|bcci;
	
	%wipe out invalid indices
	bcw = bcw(valind);
	bce = bce(valind);
	bcs = bcs(valind);
	bcn = bcn(valind);
	bcc = bcc(valind);
end

function mat = spdiag(v,k)
	
	if(nargin == 1 || k == 0)
		vsz = numel(v);
		i = (1:vsz)';
		mat = sparse(i,i,v,vsz,vsz);
	else
		vsz = numel(v);
		absk = abs(k);

		i = (1:vsz)';
		j = (1:vsz)';

		if(k>0)
			j = j+k;
		else
			i = i+absk;
		end

		mat = sparse(i,j,v,vsz+absk,vsz+absk);
	end
	
end

