function [ maxd ] = glydat(file,radius,S,flag)
% Glycosidic linkage density --> basic statistic from CircStat
% 
% INPUT VARIABLES:
%       file  - is the rescaled linkage data whose
%               1st	column - time
%               2nd ~    column - phi,psi,omega
%               last     column - free energy
%               * name without its directory
%       r     - radius of calculation area, r * r square box centered at (x,y)
%       S     - is the number of seeds used to generate the plot
%       flag  - file format options:
%               1 - plot w/o omega
%               0 - plot w/  omega
% Copy-left: Cindy Lee,2014

%% Data loading and preprocessing
clc
path(path,'~/Program/MATLAB/CircStat2012a')
str    = regexp(file,'\.','split'); name   = char(str(1));				% Disaccharide Name
dir    = strcat('~/Data/dimer/',name,'/');								% File Directory	
lch    = [dir,name,'.',num2str(S),'.linkCircHist.eps'];					% Linkage Circular Histogram
dm     = [dir,name,'.',num2str(S),'.dmatrix.',num2str(radius),'.dat'];	% Density Matrix
M      = load(strcat(dir,file));										% Data set loading
N	   = 36;															% Histogram Bin number
W      = size(M,2);	  L = S * 1e05;										% Data size
%t     = M(1:L,1);    z = M(1:L,W);										% Time and energy
x      = M(1:L,2);    y = M(1:L,3);										% Phi and Psi
x_rad  = circ_ang2rad(x); y_rad = circ_ang2rad(y);
%% Specific Data Sets Parameters 
if flag
	%% Visualization
	figure(1)
	subplot(1,2,1);
	circ_plot(x_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\phi circular distribution');	hold on;
	subplot(1,2,2);
	circ_plot(y_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\psi circular distribution');
	print(lch,'-depsc2');	hold off;	close all;
	%% Data Collection 
	dd		= mydensity([x y],'circular',radius);
	maxd	= max(dd);	nd = dd/maxd;
	data1	= [M(1:L,1:3) nd M(1:L,W)];
	dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
else 
	%% Initialization
	w	= M(1 : L, W - 1);	w_rad = circ_ang2rad(w);
	%% Visualization
	figure(1)
	subplot(1,3,1);
	circ_plot(x_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\phi circular distribution');	hold on;
	subplot(1,3,2);
	circ_plot(y_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\psi circular distribution');	hold on;
	subplot(1,3,3);
	circ_plot(w_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\omega circular distribution');
	print(lch,'-depsc2');	hold off;	close all;
	%% Data Collection 
	dd		= mydensity([x y w],'circular3',radius);		% Density measured in volume
	maxd	= max(dd);	nd	= dd/maxd;
	dd1		= mydensity([x y],'circular',radius);
	maxd1	= max(dd1)	
    nd1	= dd1/max(dd1);
	dd2		= mydensity([y w],'circular',radius);
	maxd2	= max(dd2) 	
    nd2	= dd2/max(dd2);
	data1	= [M(1:L,1:4) nd1 nd2 nd M(1:L,W)];
	dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
end %if
end %function
