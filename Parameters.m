function par = Parameters
	% Parameters
	%
	% Definite parameters go here
	
	addpath('control');
	addpath('solvers');
	addpath('bcfunctions');
	addpath('rhfunctions');
	addpath('maps');
	
	par.maptype = 'g';
	par.mapfile = 'symch.txt';
	par.h = 0.1;
	par.toPlot = 2;			%1==surf,2==quiver,3==scatter
	par.filter = true;
	par.numfilter = 1;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCSymCh; 
	par.solver = @SODuo;
	
end

