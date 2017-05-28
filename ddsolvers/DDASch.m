function psimesh = DDASch(grids,bcinds,rhs,filterMat,h)
	%DDASCH additive schwarz
	
	psimesh = SODuo(numel(grids{1}),numel(grids{2}),bcinds,rhs,filterMat,h);
	
	%TODO complete the bcinds
	[newgrids,newpsis,newbcinds] = Decompose(grids,psimesh,par.ddbounds);
	
	psimesh1 = newpsis{1};
	psimesh2 = newpsis{2};
	psimesh3 = newpsis{3};
	
	
end

