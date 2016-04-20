function [ maxdn,maxdr] = glydat3(file,r,S)
% Glycosidic linkage density --> basic statistic from CircStat
% 
% INPUT VARIABLES:
%       file  - is the rescaled linkage data whose
%               1st	column - time
%               2nd ~    column - phi,psi,omega
%               last     column - free energy
%               * name without its directory
%       r	  - radius of calculation area, r * r square box centered at (x,y)
%       S     - is the number of seeds used to generate the plot 
% Copy-left: Cindy Lee,2014

%% Data loading and preprocessing
clc
path(path,'~/Program/MATLAB/CircStat2012a')
str	   = regexp(file,'\.','split');	name = char(str(1));				% Trisaccharide Name
dir    = strcat('~/Data/trimer/',name,'/');								% File Directory	
link1  = str2double(name(1)) + 10;										% Linkage Info
link2  = str2double(name(2)) + 10;									
lch1   = [dir,name,'.',num2str(S),'.linkCircHist.n.eps'];    			% Cir Hist for non-reducing end
lch2   = [dir,name,'.',num2str(S),'.linkCircHist.r.eps'];    			% Cir Hist for reducing end
dm     = [dir,name,'.',num2str(S),'.dmatrix.',num2str(r),'.dat'];		% Density Matrix

N	   = 36;															% Histogram Bin number
M      = load(strcat(dir,file));										% Data set loading
W      = size(M,2);	  L    = S * 1e05;									% Data size
%t     = M(1:L,1);    z    = M(1:L,W);									% Time and energy
x_n	   = M(1:L,2);	xn_rad = circ_ang2rad(x_n);		                    % Phi of link 1
y_n    = M(1:L,3);	yn_rad = circ_ang2rad(y_n);                         % Psi of link 1
if link1 == 16
	z_n = M(1:L,4);	zn_rad = circ_ang2rad(z_n);							% Omg of link 1(possible)
end % link1

if link2 == 16
	x_r = M(1:L,W-3); xr_rad = circ_ang2rad(x_r);						% Phi of link 2 
	y_r = M(1:L,W-2); yr_rad = circ_ang2rad(y_r);						% Psi of link 2
	z_r = M(1:L,W-1); zr_rad = circ_ang2rad(z_r);						% Omg of link 2
else
	x_r = M(1:L,W-2); xr_rad = circ_ang2rad(x_r);						% Phi of link 2 
	y_r = M(1:L,W-1); yr_rad = circ_ang2rad(y_r);						% Psi of link 2
end % link2

%% Visualization
if link1 == 16
	figure(1)
	subplot(1,3,1);
	circ_plot(xn_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\phi circular distribution');	hold on;
	subplot(1,3,2);
	circ_plot(yn_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\psi circular distribution');	hold on;
	subplot(1,3,3);
	circ_plot(zn_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\omega circular distribution');
	print(lch1,'-depsc2');	hold off;	close all;
else
	figure(1)
    subplot(1,2,1);
    circ_plot(xn_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
    title('\phi circular distribution');    hold on;
    subplot(1,2,2);
    circ_plot(yn_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
    title('\psi circular distribution'); suptitle('circular distribution for non-reducing sugar');
    print(lch1,'-depsc2');  hold off;   close all;
end % End of link1
if link2 == 16
	figure(2)
	subplot(1,3,1);
	circ_plot(xr_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\phi circular distribution');	hold on;
	subplot(1,3,2);
	circ_plot(yr_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\psi circular distribution');	hold on;
	subplot(1,3,3);
	circ_plot(zr_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\omega circular distribution');
	print(lch2,'-depsc2');	hold off;	close all;
else
	figure(2)
	subplot(1,2,1);
	circ_plot(xr_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\phi circular distribution');	hold on;
	subplot(1,2,2);
	circ_plot(yr_rad,'hist',[],N,true,true,'linewidth',1,'color',[1 0.4 0.6]);
	title('\psi circular distribution'); suptitle('circular distribution for reducing sugar');
	print(lch2,'-depsc2');	hold off;	close all;
end % End of link2
%% Data Collection 
if link1 == 16 && link2 == 16
	dn_1	= mydensity([x_n y_n],'circular',r);	maxdn_1 = max(dn_1);	ndn_1 = dn_1/maxdn_1;
	dn_2	= mydensity([y_n z_n],'circular',r);	maxdn_2 = max(dn_2);	ndn_2 = dn_2/maxdn_2;
	dr_1	= mydensity([x_r y_r],'circular',r);	maxdr_1 = max(dr_1);	ndr_1 = dr_1/maxdr_1;
	dr_2	= mydensity([y_r z_r],'circular',r);	maxdr_2 = max(dr_2);	ndr_2 = dr_2/maxdr_2;
	maxdn	= [maxdn_1 maxdn_2];	maxdr = [maxdr_1 maxdr_2];
	ndn		= [ndn_1 ndn_2];	ndr = [ndr_1 ndr_2];
elseif link1 ~= 16 && link2 == 16
	ddn    = mydensity([x_n y_n],'circular',r); maxdn  = max(ddn);  ndn = ddn/maxdn;
	dr_1 = mydensity([x_r y_r],'circular',r);   maxdr_1 = max(dr_1);    ndr_1 = dr_1/maxdr_1;
    dr_2 = mydensity([y_r z_r],'circular',r);   maxdr_2 = max(dr_2);    ndr_2 = dr_2/maxdr_2;
	maxdr = [maxdr_1 maxdr_2];	ndr = [ndr_1 ndr_2];
elseif link1 == 16 && link2 ~= 16 
	dn_1 = mydensity([x_n y_n],'circular',r);   maxdn_1 = max(dn_1);    ndn_1 = dn_1/maxdn_1;
    dn_2 = mydensity([y_n z_n],'circular',r);   maxdn_2 = max(dn_2);    ndn_2 = dn_2/maxdn_2;
	ddr    = mydensity([x_r y_r],'circular',r); maxdr  = max(ddr);  ndr = ddr/maxdr;
	maxdn = [maxdn_1 maxdn_2];	ndn = [ndn_1 ndn_2];
else
	ddn    = mydensity([x_n y_n],'circular',r);	maxdn  = max(ddn);	ndn = ddn/maxdn;
	ddr    = mydensity([x_r y_r],'circular',r); maxdr  = max(ddr);	ndr = ddr/maxdr;
end % End of link types
	data1  = [M(1:L,1:W -1) ndn ndr M(1:L,W)];
	dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
end %function
