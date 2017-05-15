function [fx,fy,inds] = GetAppliedForce(xinit,yinit,xmesh,ymesh,valInd,sz)
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
	inds = fy;
	
	for i=1:numel(yinit)
		ind1 = xmesh==zeros(sz,1) & ymesh==yinit(i);
		ind = ind1 & valInd(ind1);
		fx(ind) = 10;
		inds = inds | ind;
	end
	
	for i=1:numel(yinit)
		ind1 = xmesh==max(xinit)*ones(sz,1) & ymesh==yinit(i);
		ind = ind1 & valInd(ind1);
		fx(ind) = 10;
		inds = inds | ind;
	end
	
	
end

