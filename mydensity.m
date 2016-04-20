function [dd] = mydensity(A,method,r)
% Computes the data density (points/area) of scattered points
% Attention: If data is circular instead of linear please use method 'circular'
%
% USAGE:
%	dd = datadensity(A,method,radius)
%
% INPUT:
%	A       - coordinates of points: x = A(:,1) y = A(:,2)
%           - coordinates of 3D points: x = A(:,1) y = A(:,2) z = A(:,3)
%	method  - either 'squares','circles', or 'voronoi','circular','circular3'
%		    - default = 'voronoi'
%	r       - Equal to the circle radius or half the square width
% Copy-left: Cindy Lee,2014
% Data preprocess
if size(A,2) == 2
	method = 'circular';
	x = A(:,1); y = A(:,2);
else
	method = 'circular3';
	x = A(:,1); y = A(:,2); z = A(:,3);
end % if
Ld = length(x);
dd = zeros(Ld,1);
%Data Density Calculation
switch method 	
	case 'squares'  %----- Using squares -----
		for k = 1 : Ld
            	dd(k) = sum(x > (x(k) - r) & x < (x(k) + r) ...
		& y > (y(k) - r) & y < (y(k) + r) );
        end %for
        area = (2*r)^2;
        dd = dd/area;
	case 'circles'
        for k = 1 : Ld
            	dd(k) = sum( sqrt((x-x(k)).^2 + (y-y(k)).^2) < r );
        end %for
        area = pi*r^2;
        dd = dd/area;
	case 'voronoi'  %----- Using voronoi cells ------
        [v,c] = voronoin([x,y]);     
        for k = 1 : length(c) 
		%If at least one of the indices is 1, then it is an open region, 
		%its area is infinity and the data density is 0
			if all(c{k}>1)   
				a = polyarea(v(c{k},1),v(c{k},2));
				dd(k) = 1/a;
			end %if
        end %for
	case 'circular' %----- Considering data's speciality -----
		for k = 1 : Ld
        	dd(k) = sum(cos((x - x(k))/180 * pi) > cos(r/360 * pi) ...
			& cos((y - y(k))/180 * pi) > cos(r/360 * pi) );
		end %for
        area = r^2;
        dd = dd/area;
	case 'circular3' %----- 3D circular data -----
		for k = 1 : Ld
			dd(k) = sum(cos((x - x(k))/180 * pi) > cos(r/360 * pi) ...
			& cos((y - y(k))/180 * pi) > cos(r/360 * pi) ...
			& cos((z - z(k))/180 * pi) > cos(r/360 * pi) );
		end %for
		volume = r^3;
		dd = dd/volume;
	end %switch
return
end % function
