xinit = (1:100)';
yinit = (1:50)';
xsz = numel(xinit);
ysz = numel(yinit);
xmesh = kron(ones(ysz,1),xinit);
ymesh = kron(yinit,ones(xsz,1));
Xmesh = reshape(xmesh,[ysz,xsz]);
Ymesh = flipud((reshape(ymesh,[ysz,xsz])));
pmesh = xmesh+ymesh;
scatter3(xmesh,ymesh,pmesh)
scatter3(xmesh,ymesh,pmesh,[],'.')
Pmesh = Xmesh+Ymesh;
figure(2)
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none','facecolor','interp')
Ymesh = reshape(ymesh,[ysz,xsz]);
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none','facecolor','interp')
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none')
Pmesh = Xmesh+Ymesh;
surf(Xmesh,Ymesh,Pmesh,'edgecolor','none')