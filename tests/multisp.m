n = 100;
a = rand;
b = rand(n);

tic
for i=1:100000
	a*b;
	if(mod(i,10000)==0)
		disp(i)
	end
end
toc

tic
for i=1:100000
	a.*b;
	if(mod(i,10000)==0)
		disp(i)
	end
end
toc

