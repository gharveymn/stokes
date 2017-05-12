function [fx,fy] = GetAppliedForce(x,y,sz)
	%APPFORCE Summary of this function goes here
	%   Detailed explanation goes here
	
	par = Parameters;
	
	file = fopen(par.initflowfile, 'r');
	formatSpec = '%f';
	data = fscanf(file,formatSpec);
	
	fx = zeros(sz,1);
	fy = fx;
	
	for i=1:4:numel(data)
		fx(x(:)==data(i)&y(:)==data(i+1)) = data(i+2);
		fy(x(:)==data(i)&y(:)==data(i+1)) = data(i+3);
	end
	
	
end

