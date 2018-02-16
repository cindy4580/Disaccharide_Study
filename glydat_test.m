function [ maxd ] = glydat_test(M,name,radius,S,T,start,incr,flag)
% Glycosidic linkage density --> basic statistic from CircStat
% 
% INPUT VARIABLES:
%       M     - Specific Data set
%       name  - Disaccharide Name (without its directory)
%       radius- Radius of area of interest
%       S     - Random number generation seed 
%       T     - Simulation temperature
%       start - Equilibrium cutoff
%       incr  - Resampling step size
%       flag  - file format options:
%               1 - plot w/o omega
%               0 - plot w/  omega
% Note: Input data in degrees [0 360] instead of radians
% Copy-left: Xindi Li, 2014
 
dir = strcat('~/Data/dimer/',name,'/');					   
dm = [dir,'DCD/MC2long/',name,'.',num2str(S),'.',num2str(T),...
      '.start',num2str(start),'.skip',num2str(incr),'.d.dat'];
    

W  = size(M,2);	
x  = M(1:incr:end,2);    y = M(1:incr:end,3);							  

if flag
	dd		= mydensity([x y],'circular',radius);
	nd 		= dd/max(dd);
	data1	= [M(1:incr:end, 1:3) nd M(1:incr:end, W)];
	dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
else 
	w		= M(1:incr:end, W-1);
	dd		= mydensity([x y w],'circular3',radius);	   % Volume density
	maxd	= max(dd);	nd	= dd/maxd;
	dd1		= mydensity([x y],'circular',radius);
	maxd1	= max(dd1)	
    nd1	= dd1/max(dd1);
	dd2		= mydensity([y w],'circular',radius);
	maxd2	= max(dd2) 	
    nd2	= dd2/max(dd2);
	data1	= [M(1:incr:end,1:4) nd1 nd2 nd M(1:incr:end,W)];
	dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
end %if
end %function
