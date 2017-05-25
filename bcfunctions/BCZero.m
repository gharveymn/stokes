function [rhs,bcinds] = BCZero(xmesh,ymesh,rhs,on)
	bcinds = on;
	rhs(bcinds) = 0;
end