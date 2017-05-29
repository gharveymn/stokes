xinit = (1:100)';
yinit = (1:50)';
nx = numel(xinit);
ny = numel(yinit);
xmesh = kron(ones(ny,1),xinit);
ymesh = kron(yinit,ones(nx,1));
Xmesh = reshape(xmesh,[ny,nx]);
Ymesh = flipud((reshape(ymesh,[ny,nx])));
pmesh = xmesh+ymesh;
scatter3(xmesh,ymesh,pmesh)
scatter3(xmesh,ymesh,pmesh,[],'.')
Pmesh = Xmesh+Ymesh;
figure(2)
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none','facecolor','interp')
Ymesh = reshape(ymesh,[ny,nx]);
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none','facecolor','interp')
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none')
Pmesh = Xmesh+Ymesh;
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none')