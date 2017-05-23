function Plot(mat,vec,toPlot)
	
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
	
	if(toPlot == "surf")
		X = X(2:end-1,2:end-1);
		Y = Y(2:end-1,2:end-1);
		U = U(2:end-1,2:end-1);
		V = V(2:end-1,2:end-1);
		Psi = Psi(2:end-1,2:end-1);
		
		figure(1)
		surf(X,Y,U,'edgecolor','none','facecolor','interp')
		axis(ax)

		figure(2)
		surf(X,Y,V,'edgecolor','none','facecolor','interp')
		axis(ax)

		figure(3)
		surf(X,Y,Psi,'edgecolor','none','facecolor','interp')
		axis(ax)
	elseif(toPlot == "quiver")
		figure(1)
		ax = MakeAxis(Xmesh,Ymesh);
		quiver(xmesh,ymesh,u,v)
		axis(ax)

		figure(3)
		ax = MakeAxis(Xmesh,Ymesh);
		clrs = abs(psi)./max(abs(psi(isfinite(psi))));
		scatter3(xmesh,ymesh,psi,10,[clrs zeros(numel(clrs),1) 1-clrs], '.')
		axis(ax)
	else
		figure(1)
		scatter3(x,y,u,10,'.')
		axis(ax)

		figure(2)
		scatter3(x,y,v,10,'.')
		axis(ax)

		figure(3)
		scatter3(x,y,psi,10,'.')
		axis(ax)
	end
	
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



