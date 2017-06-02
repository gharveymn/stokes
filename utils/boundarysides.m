function [bc,bcfull] = boundarysides(xmeshfull,ymeshfull,onfull,valind,nx)
	%GETWHEREBOUNDARIES I'm somewhat suprised this actually works
	
	xmin = min(xmeshfull);
	xmax = max(xmeshfull);
	ymin = min(ymeshfull);
	ymax = max(ymeshfull);
	
	%set to contain all possible indices
	bcw = onfull;
	bce = onfull;
	bcs = onfull;
	bcn = onfull;
	
	%get indices on the ends
	xminb = (xmeshfull==xmin);
	xmaxb = (xmeshfull==xmax);
	yminb = (ymeshfull==ymin);
	ymaxb = (ymeshfull==ymax);
	
	%if an index is on the overall boundary then it must be a boundary point
	bcw = bcw&xminb;
	bce = bce&xmaxb;
	bcs = bcs&yminb;
	bcn = bcn&ymaxb;
	
	%make shifted index matrices, filter out those indices at the max values since circshift loops
	r = circshift(valind&~xmaxb,1);
	l = circshift(valind&~xminb,-1);
	u = circshift(valind&~ymaxb,nx);
	d = circshift(valind&~yminb,-nx);
	
	%if the shift makes it go off the boundary then we have a direction
	bcw = bcw|(onfull&~r);
	bce = bce|(onfull&~l);
	bcs = bcs|(onfull&~u);
	bcn = bcn|(onfull&~d);
	
	%corners--is in two of the previous or is surrounded
	bcc = (bcw&bce)|(bcw&bcs)|(bcw&bcn)|(bce&bcs)|(bce&bcn)|(bcs&bcn);
	
	%inner corner boundary condition
	bcci = onfull&(r&l&u&d);
	
	%bcw = bcw|bcci;
	%bce = bce|bcci;
	%bcs = bcs|bcci;
	%bcn = bcn|bcci;
	bcc = bcc|bcci;
	
	bcfull = {bcw,bce,bcs,bcn,bcc};
	
	%wipe out invalid indices
	bcw = bcw(valind);
	bce = bce(valind);
	bcs = bcs(valind);
	bcn = bcn(valind);
	bcc = bcc(valind);
	
	bc = {bcw,bce,bcs,bcn,bcc};
end

