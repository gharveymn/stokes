classdef Solver
	%SOLVER Provides an interface for all solvers
	
	properties
		grids
		filtering
		bcinds
		rhs
		par
		mats
	end
	
	methods (abstract)
		solve(grids{9},grids{10},bcinds,rhs,filtering{1},par.h,mats)
	end
	
end

