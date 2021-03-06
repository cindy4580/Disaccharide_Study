function [ maxd ] = glydat(name,radius,S,T,start,incr,flag)
% Glycosidic linkage density --> basic statistic from CircStat
% 
% INPUT VARIABLES:
%       name  - Disaccharide Name (without its directory)
%       r     - Radius of area of interest
%       S     - Random number generation seed 
%       T     - Simulation temperature
%       start - Equilibrium cutoff
%       incr  - Resampling step size
%       flag  - file format options:
%               1 - plot w/o omega
%               0 - plot w/  omega
% Note: Input data in degrees [0 360] instead of radians
% Copy-left: Xindi Li, 2014

clc; 
dir = strcat('~/Data/dimer/',name,'/');					   
% M      = load(strcat(dir,'orig_data/',name,'.',num2str(S),...
%               '.',num2str(T),'.pl.dat'));
% dm     = [dir,'density/',name,'.',num2str(S),'.',...
%           num2str(T),'.d',num2str(incr),'.dat']; 
M = load(strcat(dir,'DCD/DMAX30/',name,'.',num2str(S),...
    '.',num2str(T),'.pl.dat'));
dm = [dir,'DCD/DMAX30/',name,'.',num2str(S),'.',num2str(T),...
    '.d',num2str(incr),'.dat'];
    

W  = size(M,2);	
x  = M(start:incr:end,2);    y = M(start:incr:end,3);							  

if flag
	dd		= mydensity([x y],'circular',radius);
	nd 		= dd/max(dd);
	data1	= [M(start:incr:end, 1:3) nd M(start:incr:end, W)];
	dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
else 
	w		= M(start:incr:end, W-1);
	dd		= mydensity([x y w],'circular3',radius);	   % Volume density
	maxd	= max(dd);	nd	= dd/maxd;
	dd1		= mydensity([x y],'circular',radius);
	maxd1	= max(dd1)	
    nd1	= dd1/max(dd1);
	dd2		= mydensity([y w],'circular',radius);
	maxd2	= max(dd2) 	
    nd2	= dd2/max(dd2);
	data1	= [M(start:incr:end,1:4) nd1 nd2 nd M(start:incr:end,W)];
	dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
end %if
end %function
