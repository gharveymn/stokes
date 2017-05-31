function par = Parameters
	% Parameters
	%
	% Definite parameters go here
	
	addpath('ddsolvers');
	addpath('solvers');
	addpath('bcfunctions');
	addpath('rhfunctions');
	addpath('maps');
	addpath('utils');
	
	par.maptype = 'g';
	par.mapfile = 'symch.txt';
	par.h = 0.05;
	par.toPlot = 4;			%1==surf,2==quiver,3==scatter,4==contour
	par.filter = true;
	par.numfilter = 1;
	par.ghostpoints = true;
	par.zeroout = false;
	
	%domain decomposition parameters
	par.ddrun = false;
	par.ddbounds = {{[0.0,-0.5],[1.5,0.5]},{[1.0,-1.5],[3.5,1.5]},{[3.0,-1.5],[5.0,1.5]}};
	par.ddoverlap = 0.5;
	par.ddmidratio = 0.6;
	par.dditer = 5;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCSymCh; 
	par.solver = @SODuo;
	par.ddsolver = @DDMSch;
	
end

