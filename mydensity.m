function [dd] = mydensity(A,method,r)
% Computes the data density (points/area) of scattered points
% Note: If data is circular instead of linear, use method 'circular';
%       Data sets are in degrees instead of radians;
%
% USAGE:
%	dd = datadensity(A,method,radius)
%
% INPUT:
%	A       - coordinates of points: x = A(:,1) y = A(:,2)
%           - coordinates of 3D points: x = A(:,1) y = A(:,2) z = A(:,3)
%	method  - 'squares','circles','voronoi','circular',or 'circular3'
%	r       - Equal to the circle radius or half the square width
% Copy-left: Cindy Lee,2014

%% Data preprocess
if strcmp(method,'circular') && size(A,2) == 3
   method = 'circular3';
end
x  = A(:,1);    y = A(:,2);
Ld = length(x);
dd = zeros(Ld,1);
%% Data Density Calculation
switch method 	
	case 'squares'                          %----- Using square boxes -----
        for k = 1 : Ld
            dd(k) = sum(x > (x(k) - r) & x < (x(k) + r) ...
                      & y > (y(k) - r) & y < (y(k) + r) );
        end
        a   = (2*r)^2;
        dd  = dd/a;
	case 'circles'
        for k = 1 : Ld
            dd(k) = sum( sqrt((x-x(k)).^2 + (y-y(k)).^2) < r );
        end 
        a   = pi*r^2;
        dd  = dd/a;
	case 'voronoi'                        %----- Using voronoi cells ------
        [v,c] = voronoin([x,y]);     
        for k = 1 : length(c) 
		%If at least one of the indices is 1, then it is an open region, 
		%its area is infinity and the data density is 0
			if all(c{k}>1)   
				a = polyarea(v(c{k},1),v(c{k},2));
				dd(k) = 1/a;
			end %if
        end %for
	case 'circular'                          %----- Circular Data: 2D -----
		for k = 1 : Ld
        	dd(k) = sum(cos((x - x(k))/180 * pi) > cos(r/360 * pi) ...
                      & cos((y - y(k))/180 * pi) > cos(r/360 * pi) );
		end %for
        a   = r^2;
        dd  = dd/a;
	case 'circular3'                         %----- Circular data: 3D -----
		z = A(:,3);
        for k = 1 : Ld
			dd(k) = sum(cos((x - x(k))/180 * pi) > cos(r/360 * pi) ...
                      & cos((y - y(k))/180 * pi) > cos(r/360 * pi) ...
                      & cos((z - z(k))/180 * pi) > cos(r/360 * pi) );
		end %for
		vl = r^3;
		dd = dd/vl;
	end %switch
return
end % function
