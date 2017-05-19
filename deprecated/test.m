minx = -1;
maxx = 1;
miny = -2;
maxy = 2;
h = 0.01;

xinit = (minx:h:maxx)';
yinit = (miny:h:maxy)';

xsz = numel(xinit);
ysz = numel(yinit);

[X,Y] = meshgrid(xinit,yinit);
X = X';
Y = Y';

L_h = LaplacianFactory(xsz,ysz,h);

x = kron(ones(ysz,1),xinit);
y = kron(yinit,ones(xsz,1));
fxy = x.^4 + y.^3;
d2fxy = L_h*fxy;
d2fXY = reshape(d2fxy,numel(xinit),numel(yinit));

syms xi yi
f(xi,yi) = xi^4 + yi^3;
df(xi,yi) = diff(f(xi,yi),xi,2) + diff(f(xi,yi),yi,2);

figure(1)
fsurf(df,[minx,maxx,miny,maxy],'edgecolor','none','facecolor','interp')

figure(2)
surf(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),d2fXY(2:end-1,2:end-1),'edgecolor','none','facecolor','interp')
%scatter3(x,y,d2fxy,[],'.')

