function [bcx,bcy,inds] = GetDirichBC(xinit,yinit,xmesh,ymesh,valInd,sz)
	
	bcx = zeros(sz,1);
	bcy = bcx;
	inds = bcy;
	
	for i=1:numel(yinit)
		ind1 = xmesh==zeros(sz,1) & ymesh==yinit(i);
		ind = ind1 & valInd(ind1);
		bcx(ind) = 4*0.3*yinit(i)*(0.5 - yinit(i))/(0.41)^2;
		inds = inds | ind;
	end
	
% 	for i=1:numel(yinit)
% 		ind1 = xmesh==max(xinit)*ones(sz,1) & ymesh==yinit(i);
% 		ind = ind1 & valInd(ind1);
% 		fx(ind) = 10;
% 		inds = inds | ind;
% 	end
	
	
end

