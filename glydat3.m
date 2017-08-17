function [ ] = glydat3(name,r,S,varargin)
% Glycosidic linkage density --> basic statistic from CircStat
% 
% INPUT VARIABLES:
%       name  - Trisaccharide Name (without its directory)
%       r	  - Radius of area centered at (x,y)
%       S     - Random number generation seed 
%
% Copy-left: Cindy Lee,2014

%% Data loading and preprocessing
clc;
path(path,'~/Program/MATLAB/CircStat2012a')
dir     = strcat('~/Data/trimer/',name,'/');				%File Directory
l1      = str2double(name(1)) + 10;							 % Linkage Info
l2      = str2double(name(2)) + 10;									
% lch1    = [dir,name,'.',num2str(S),'.linkCircHist.n.eps']; 
% lch2    = [dir,name,'.',num2str(S),'.linkCircHist.r.eps'];   
dm      = [dir,'density/',name,'.',num2str(S),'.d.dat'];		
N       = 36;										 % Histogram Bin number
M       = load(strcat(dir,'orig_data/',name,'.',num2str(S),'.pl.dat'));
W       = size(M,2);     L   = 1;								% Data size

% Non-reducing glycosidic linkage
x_n     = M(1:L:end,2);  y_n = M(1:L:end,3);
if l1   == 16
	z_n = M(1:L:end,4);							
end
% Reducing glycosidic linkage
if l2   == 16
	x_r = M(1:L:end,W-3);   y_r = M(1:L:end,W-2);   z_r = M(1:L:end,W-1);
else
	x_r = M(1:L:end,W-2);   y_r = M(1:L:end,W-1);						 
end

%% Circular Distribution: Non-reducing End
% if l1 == 16
% 	figure(1)
% 	subplot(1,3,1);
% 	circ_plot(circ_ang2rad(x_n),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
% 	title('\phi circular distribution');	hold on;
% 	subplot(1,3,2);
% 	circ_plot(circ_ang2rad(y_n),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
% 	title('\psi circular distribution');	hold on;
% 	subplot(1,3,3);
% 	circ_plot(circ_ang2rad(z_n),'hist',[],N,true,true,'linewidth',1,...
%         'color',[1 0.4 0.6]);
% 	title('\omega circular distribution');
%     suptitle('Circular Distribution for Non-reducing sugar');
% 	print(lch1,'-depsc2');	hold off;	close all;
% else
% 	figure(1)
%     subplot(1,2,1);
%     circ_plot(circ_ang2rad(x_n),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
%     title('\phi circular distribution');    hold on;
%     subplot(1,2,2);
%     circ_plot(circ_ang2rad(y_n),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
%     title('\psi circular distribution'); 
%     suptitle('Circular Distribution for Non-reducing Sugar');
%     print(lch1,'-depsc2');  hold off;   close all;
% end 
%% Circular Distribution: Reducing End
% if l2 == 16
% 	figure(2)
% 	subplot(1,3,1);
% 	circ_plot(circ_ang2rad(x_r),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
% 	title('\phi circular distribution');	hold on;
% 	subplot(1,3,2);
% 	circ_plot(circ_ang2rad(y_r),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
% 	title('\psi circular distribution');	hold on;
% 	subplot(1,3,3);
% 	circ_plot(circ_ang2rad(z_r),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
% 	title('\omega circular distribution');
%     suptitle('Circular Distribution for Reducing Sugar');
% 	print(lch2,'-depsc2');	hold off;	close all;
% else
% 	figure(2)
% 	subplot(1,2,1);
% 	circ_plot(circ_ang2rad(x_r),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
% 	title('\phi circular distribution');	hold on;
% 	subplot(1,2,2);
% 	circ_plot(circ_ang2rad(y_r),'hist',[],N,true,true,'linewidth',1,...
%              'color',[1 0.4 0.6]);
% 	title('\psi circular distribution'); 
%     suptitle('Circular Distribution for Reducing Sugar');
% 	print(lch2,'-depsc2');	hold off;	close all;
% end
%% Data Collection: density matrix
if l1 == 16 && l2 == 16
	dn_1	= mydensity([x_n y_n],'circular',r);    
	dn_2	= mydensity([y_n z_n],'circular',r);
	dr_1	= mydensity([x_r y_r],'circular',r);	 
	dr_2	= mydensity([y_r z_r],'circular',r);	 
	ndn		= [dn_1/max(dn_1) dn_2/max(dn_2)];	
    ndr     = [dr_1/max(dr_1) dr_2/max(dr_2)];
elseif l1  ~= 16 && l2 == 16
	ddn     = mydensity([x_n y_n],'circular',r);    ndn   = ddn/max(ddn);
	dr_1    = mydensity([x_r y_r],'circular',r);     
    dr_2    = mydensity([y_r z_r],'circular',r);     
    ndr     = [dr_1/max(dr_1) dr_2/max(dr_2)];
elseif l1   == 16 && l2 ~= 16 
	dn_1    = mydensity([x_n y_n],'circular',r);     
    dn_2    = mydensity([y_n z_n],'circular',r);     
	ddr     = mydensity([x_r y_r],'circular',r);    ndr   = ddr/max(ddr);
    ndn     = [dn_1/max(dn_1) dn_2/max(dn_2)];
else
	ddn     = mydensity([x_n y_n],'circular',r);     ndn = ddn/max(ddn);
	ddr     = mydensity([x_r y_r],'circular',r);     ndr = ddr/max(ddr);
end % End of link types

data1       = [M(1:L:end,1:W -1) ndn ndr M(1:L:end,W)];
dlmwrite(dm,data1,'delimiter','\t','newline','Unix');
end %function
