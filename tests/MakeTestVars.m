h = 0.01;
xinit = (0:h:1)';
yinit = xinit;

nx = numel(xinit);
ny = numel(yinit);

xmesh = kron(ones(ny,1),xinit);
ymesh = kron(yinit,ones(nx,1));
Xmesh = (reshape(xmesh,[nx,ny]))';
Ymesh = (reshape(ymesh,[ny,nx]))';