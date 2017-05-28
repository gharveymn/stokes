function [rhs,bcinds] = BCZero(xmesh,ymesh,rhs,on,del,par)
	bcinds = on;
	rhs(bcinds) = 0;
end