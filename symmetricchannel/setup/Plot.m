function figs = Plot(mat,vec,toPlot,filter,figs)
	
	if(nargin == 4)
		figs = InitialPlot(mat,vec,toPlot,filter);
	else
		%update if we already initialized the plots
		Update(mat,vec,toPlot,filter,figs)
	end
end

function figs=InitialPlot(mat,vec,toPlot,filter)
	
	X = mat(:,:,1);
	Y = mat(:,:,2);
	U = mat(:,:,3);
	V = mat(:,:,4);
	Psi = mat(:,:,5);
	
	x = vec(:,1);
	y = vec(:,2);
	u = vec(:,3);
	v = vec(:,4);
	psi = vec(:,5);
	
	ax = MakeAxis(x,y);
	
	if(filter)
		par = Parameters;
		nf = par.numfilter;
		
		[isz,jsz] = size(X);
		z = zeros(nf,1);
		ifil = kron(ones(jsz,1),[z;ones(isz-2*nf,1)-2;z]);
		jfil = kron([z;ones(jsz-2*nf,1);z],ones(isz,1));
		fil = ifil&jfil;
		x = x(fil);
		y = y(fil);
		u = u(fil);
		v = v(fil);
		psi = psi(fil);
		
		
		X = X(1+nf:end-nf,1+nf:end-nf);
		Y = Y(1+nf:end-nf,1+nf:end-nf);
		U = U(1+nf:end-nf,1+nf:end-nf);
		V = V(1+nf:end-nf,1+nf:end-nf);
		Psi = Psi(1+nf:end-nf,1+nf:end-nf);
	end
	
	if(toPlot == "surf")
		
		figure(1)
		f1 = surf(X,Y,U,'edgecolor','none','facecolor','interp');
		axis(ax)

		figure(2)
		f2 = surf(X,Y,V,'edgecolor','none','facecolor','interp');
		axis(ax)
		
		figure(3)
		f3 = surf(X,Y,Psi,'edgecolor','none','facecolor','interp');
		axis(ax)
		
		figs = {f1,f2,f3};
		
	elseif(toPlot == "quiver")
		figure(1)
		ax = MakeAxis(X,Y);
		f1 = quiver(x,y,u,v);
		axis(ax)

		figure(2)
		ax = MakeAxis(X,Y);
		clrs = abs(psi)./max(abs(psi(isfinite(psi))));
		f2 = scatter3(x,y,psi,10,[clrs zeros(numel(clrs),1) 1-clrs], '.');
		axis(ax)
		
		figs = {f1,f2};
		
	else
		
		figure(1)
		f1 = scatter3(x,y,u,10,'.');
		axis(ax)

		figure(2)
		f2 = scatter3(x,y,v,10,'.');
		axis(ax)

		figure(3)
		f3 = scatter3(x,y,psi,10,'.');
		axis(ax)
		
		figs = {f1,f2,f3};
	end
	
	drawnow;
	
end

function Update(mat,vec,toPlot,filter,figs)
	
	U = mat(:,:,3);
	V = mat(:,:,4);
	Psi = mat(:,:,5);
	
	u = vec(:,3);
	v = vec(:,4);
	psi = vec(:,5);
	
	if(filter)
		[isz,jsz] = size(U);
		ifil = kron(ones(jsz,1),[0;ones(isz-2,1)-2;0]);
		jfil = kron([0;ones(jsz-2,1);0],ones(isz,1));
		fil = ifil&jfil;
		u = u(fil);
		v = v(fil);
		psi = psi(fil);
		
		U = U(2:end-1,2:end-1);
		V = V(2:end-1,2:end-1);
		Psi = Psi(2:end-1,2:end-1);
	end
	
	if(toPlot == "surf")
		
		set(figs{1},'ZData',U);
		set(figs{2},'ZData',V);
		set(figs{3},'ZData',Psi);
		
	elseif(toPlot == "quiver")
		
		set(figs{1},'UData',u);
		set(figs{1},'VData',v);
		
		clrs = MakeClrs(psi);
		set(figs{2},'ZData',psi);
		set(figs{2},'CData',clrs);
		
	else
		clrs = MakeClrs(u);
		set(figs{1},'ZData',u);
		set(figs{1},'CData',clrs);
		
		clrs = MakeClrs(v);
		set(figs{2},'ZData',v);
		set(figs{2},'CData',clrs);
		
		clrs = MakeClrs(psi);
		set(figs{3},'ZData',psi);
		set(figs{3},'CData',clrs);
	end
	
	drawnow;
	
end

function ax = MakeAxis(x,y)
	x = x(isfinite(x));
	y = y(isfinite(y));
	
	minx = min(min(x));
	maxx = max(max(x));
	miny = min(min(y));
	maxy = max(max(y));
	absmin = min(minx,miny);
	absmax = max(maxx,maxy);
	centerx = (maxx+minx)/2;
	centery = (maxy+miny)/2;
	difa = absmax-absmin;
	ax = [centerx-difa/2, centerx+difa/2, centery-difa/2, centery+difa/2];
end

function clrs = MakeClrs(v)
	clr = abs(v)./norm(v(isfinite(v)),inf);
	clrs = [clr zeros(numel(clr),1) 1-clr];
end



