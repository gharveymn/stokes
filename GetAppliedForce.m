function [fx,fy] = GetAppliedForce(xinit,yinit,xmesh,ymesh,valInd,sz)
	%APPFORCE Summary of this function goes here
	%   Detailed explanation goes here
	
% 	par = Parameters;
% 	
% 	file = fopen(par.initflowfile, 'r');
% 	formatSpec = '%f';
% 	data = fscanf(file,formatSpec);
% 	
% 	fx = zeros(sz,1);
% 	fy = fx;
% 	
% 	for i=1:4:numel(data)
% 		fx(xmesh(:)==data(i)&ymesh(:)==data(i+1)) = data(i+2);
% 		fy(xmesh(:)==data(i)&ymesh(:)==data(i+1)) = data(i+3);
% 	end
	
	fx = zeros(sz,1);
	fy = fx;
	
	for i=1:numel(yinit)
		fx(xmesh==0&ymesh==yinit(i)&valInd(xmesh==0&ymesh==yinit(i))) = 0.5;
	end
	
	
end

