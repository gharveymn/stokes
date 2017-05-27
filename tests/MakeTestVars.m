h = 0.01;
xinit = (0:h:1)';
yinit = xinit;

xsz = numel(xinit);
ysz = numel(yinit);

xmesh = kron(ones(ysz,1),xinit);
ymesh = kron(yinit,ones(xsz,1));
Xmesh = (reshape(xmesh,[xsz,ysz]))';
Ymesh = (reshape(ymesh,[ysz,xsz]))';