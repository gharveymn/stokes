function par = Parameters
	% Parameters
	%
	% Definite parameters go here
	
	addpath('../functions');
	addpath('../maps');
	
	par.maptype = 'g';
	par.mapfile = 'Ldomain.txt';
	par.h = 0.05;
	par.toPlot = 1;			%1==surf,2==quiver,3==scatter
	par.filter = false;
	par.numfilter = 1;
	
	par.func = @Duo;
	
end

