function par = Parameters
	% Parameters
	%
	% Definite parameters go here
	
	addpath('solvers');
	addpath('bcfunctions');
	addpath('rhfunctions');
	addpath('maps');
	
	par.maptype = 'g';
	par.mapfile = 'unit.txt';
	par.h = 0.05;
	par.toPlot = 1;			%1==surf,2==quiver,3==scatter
	par.filter = true;
	par.numfilter = 1;
	
	par.rhfunc = @RHOne;
	par.bcfunc = @BCZero; 
	par.solver = @SOSpectral;
	
end

