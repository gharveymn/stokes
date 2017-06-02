function rhs = extendgp(rhs,bcfull,valindouter,gpca,nx)
	%EXTENDGP extends rhs to ghostpoints order 
	%   Detailed explanation goes here
	
	bcw = bcfull{1};
	bce = bcfull{2};
	bcs = bcfull{3};
	bcn = bcfull{4};
	bcc = bcfull{5};
	
	%set outer regions
	
	%corners (up to second order)
	rhs = makecorners(rhs,bcc,valindouter,gpca,nx);	
	
	%west
	bcw1 = circshift(bcw,-1);
	bcw2 = circshift(bcw,-2);
	rhs(bcw1(valindouter)) = rhs(bcw(valindouter));
	rhs(bcw2(valindouter)) = rhs(bcw(valindouter));
	
	%east
	bce1 = circshift(bce,1);
	bce2 = circshift(bce,2);
	rhs(bce1(valindouter)) = rhs(bce(valindouter));
	rhs(bce2(valindouter)) = rhs(bce(valindouter));
	
	%south
	bcs1 = circshift(bcs,-nx);
	bcs2 = circshift(bcs,-2*nx);
	rhs(bcs1(valindouter)) = rhs(bcs(valindouter));
	rhs(bcs2(valindouter)) = rhs(bcs(valindouter));
	
	%north
	bcn1 = circshift(bcn,nx);
	bcn2 = circshift(bcn,2*nx);
	rhs(bcn1(valindouter)) = rhs(bcn(valindouter));
	rhs(bcn2(valindouter)) = rhs(bcn(valindouter));
	
end

function rhs = makecorners(rhs,bcc,valind,gpca,nx)
	%input i for the order
	
	for k=1:numel(gpca)
		gp = gpca{k};
		for i=-k:k
			for j=-k:k
				bccsw = circshift(circshift(bcc,-i),-j*nx) & gp;
				bccswr = circshift(circshift(bccsw,j*nx),i);

				bccnw = circshift(circshift(bcc,-i),j*nx) & gp;
				bccnwr = circshift(circshift(bccnw,-j*nx),i);

				bccse = circshift(circshift(bcc,i),-j*nx) & gp;
				bccser = circshift(circshift(bccse,j*nx),-i);

				bccne = circshift(circshift(bcc,i),j*nx) & gp;
				bccner = circshift(circshift(bccne,-j*nx),-i);

				rhs(bccsw(valind)) = rhs(bccswr(valind));
				rhs(bccnw(valind)) = rhs(bccnwr(valind));
				rhs(bccse(valind)) = rhs(bccser(valind));
				rhs(bccne(valind)) = rhs(bccner(valind));
			end
		end
	end
end

