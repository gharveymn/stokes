a = rand(4000);
b = rand(4000);
c = rand(4000);
d = rand(4000);
e = rand(4000);

bigcell = {a,b,c,d,e};

numl = 1000;

tic
for i=1:numl
	unpack(bigcell);
end
toc

tic
for i=1:numl
	direct(a);
end
toc

function unpack(sell)
	f = sell{1};
	reshape(f,size(f,1)*size(f,2),1);
end

function direct(f)
	reshape(f,size(f,1)*size(f,2),1);
end