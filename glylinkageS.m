function [ ] = glylinkageS(name,S,flag,seed)
% Glycosidic linkage --> Ramachandran plot (phi,psi,omega) for
% DISACCHARIDES of Single Seed
% Scatter plot with color indicating data density
% INPUT VARIABLES:
%       name  - is the name for specific disaccharide
%               1st	column - time
%               2nd ~    column - phi,psi,omega
%               2nd last column - normalized density
%               last     column - free energy
%               * name without its directory
%       S     - is the number of seeds used to generate the plot
%       flag  - file format options:
%               1 - plot w/o omega
%               0 - plot w/  omega
%       seed  - Specify the seed (only for one-seed plot)
% Copy-left: Cindy Lee,2014
clc
dir		= strcat('~/Data/dimer/',name,'/density/');	% File Directory
M		= load(strcat(dir,name,'.',num2str(seed),'.d.dat'));
W		= size(M,2);	L = S * 1e05;							% Data size
x		= M(1:L,2);	y = M(1:L,3);								% (Phi,Psi)
d		= M(1:L,W - 1);											% Density
dot		= M((M(:,W) == min(M(:,W))),:);							% Lowest energy matrix
color	= hsv(10);
xmax	= 360;	xmin = 0;
ymax	= 360;	ymin = 0;
xbin	= (xmin + 5) :10: (xmax - 5);	
ybin	= (ymin + 5) :10: (ymax - 5);
rn		= [dir,name,'.',num2str(seed),'.rama.eps'];	% Linkage Plot
%% Specific Data Sets Parameters
if flag
	switch name
        % ---------------- 1 -> 2 linkages --------------------------------------------
        case '12AA_MM'			% 1
		T = 'Disaccharide: \alphaMannose-\alpha(1,2)-\alphaMannose';
		case '12BA_NM'			% 2
		T = 'Disaccharide: \betaGlcNAc-\beta(1,2)-\alphaMannose';
		case '12AB_FG'          % 3
		T = 'Disaccharide: \alphaFucose-\alpha(1,2)-\betaGalactose';
		case '12BA_XM'			% 4
		T = 'Disaccharide: \betaXylose-\beta(1,2)-\betaMannose';
		% ---------------- 1 -> 3 linkages --------------------------------------------
		case '13AB_FN'			% 1
		T = 'Disaccharide: \alphaFucose-\alpha(1,3)-\betaGlcNAc';
		case '13AB_GG'			% 2
		T = 'Disaccharide: \alphaGalactose-\alpha(1,3)-\betaGalactose';
		case '13AB_MM'			% 3
		T = 'Disaccharide: \alphaMannose-\alpha(1,3)-\betaMannose';
		case '13AB_MM_A'		% 4
		T = 'Disaccharide: \alphaMannose-\alpha(1,3)-\alphaMannose';
		case '13AB_YG'			% 5
		T = 'Disaccharide: \alphaGalNAc-\alpha(1,3)-\betaGalactose';
		case '13AB_YY'			% 6
		T = 'Disaccharide: \alphaGalNAc-\alpha(1,3)-\betaGalNAc';
		%------------------------------------------------------------%
		case '13BB_GN'			% 1
		T = 'Disaccharide: \betaGalactose-\beta(1,3)-\betaGlcNAc';
		case '13BB_GY'			% 2
		T = 'Disaccharide: \betaGalactose-\beta(1,3)-\betaGalNAc';
		case '13BB_GY_A'		% 3
		T = 'Disaccharide: \betaGalactose-\beta(1,3)-\alphaGalNAc';
		case '13BB_NG'			% 4
		T = 'Disaccharide: \betaGlcNAc-\beta(1,3)-\betaGalactose';
        case '13BB_NY_A'		% 5
        T = 'Disaccharide: \betaGlcNAc-\beta(1,3)-\alphaGalNAc';
		case '13BB_UG'			% 6
        T = 'Disaccharide: \betaGlcA-\beta(1,3)-\betaGalactose';
		case '13BB_YG'			% 7
		T = 'Disaccharide: \betaGalNAc-\beta(1,3)-\alphaGalactose';
		% ---------------- 1 -> 4 linkages --------------------------------------------
		case '14AA_GG'			% 1
		T = 'Disaccharide: \alphaGalactose-\alpha(1,4)-\betaGalactose';
		case '14AB_FN'			% 2
		T = 'Disaccharide: \alphaFucose-\alpha(1,4)-\betaGlcNAc';
		%----------------------------------------------------------------%	
		case '14BB_MN'			% 1
		T = 'Disaccharide: \betaMannose-\beta(1,4)-\betaGlcNAc';
		case '14BB_NN'			% 2
		T = 'Disaccharide: \betaGlcNAc-\beta(1,4)-\betaGlcNAc';
		case '14BB_NM'			% 3
		T = 'Disaccharide: \betaGlcNAc-\beta(1,4)-\betaMannose';
		case '14BB_GP'			% 4
		T = 'Disaccharide: \betaGalactose-\beta(1,4)-\betaGlucose';
		case '14BA_YG'			% 5
		T = 'Disaccharide: \betaGalNAc-\beta(1,4)-\betaGlactose';
		case '14BB_YN'			% 6
		T = 'Disaccharide: \betaGalNAc-\beta(1,4)-\betaGlcNAc';
		case '14BA_GY'			% 7
		T = 'Disaccharide: \betaGalactose-\beta(1,4)-\betaGalNAc';
		case '14BB_GN'			% 8
		T = 'Disaccharide: \betaGalactose-\beta(1,4)-\betaGlcNAc';		
		case '14BB_NM_A'		% 9
		T = 'Disaccharide: \betaGlcNAc-\beta(1,4)-\alphaMannose';
		% ----------------------- Sialic linkages -------------------------------------
		case '23AB_SG'			% 1
		T = 'Disaccharide: \alphaNeu5Ac-\alpha(2,3)-\betaGalactose';
	end %switch
    %% Visualization
	h12  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
	nx  = hist(x,xbin);
	bar(xbin,nx/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
	xlim([xmin xmax]);
	ylabel('Frequency','FontSize',8);	set(gca,'fontsize',7);	hold on;

	h13  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
	ny  = hist(y,ybin);
	bar(ybin,ny/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
	xlim([ymin,ymax]);	view(h13,270,90);	set(gca,'xticklabel',[],'fontsize',7);
	ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);	hold on;

	h14  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
	fastscatter(x,y,d);	colorbar;	xlim([xmin xmax]); ylim([ymin ymax]);	hold on;
	for i = 1 : size(dot,1)
		plot([dot(i,2) dot(i,2)], [ymin ymax],'--','color',color(i,:));	hold on;
		plot([xmin xmax], [dot(i,3) dot(i,3)],'--','color',color(i,:));	hold on;
        plot(dot(i,2),dot(i,3),'+','color',color(i,:));
	end % for
	grid on;	box on;
	xlabel('\phi(degree)','Fontsize',12);	suptitle(T);    % No other input arguments allowed
	hold off;	print(rn,'-depsc2');
else
	switch name
		case '16AB_FN'			% 1
			T  = 'Disaccharide: \alphaFucose-\alpha(1,6)-\betaGlcNAc';
			T1 = 'Disaccharide: \alphaFucose-\alpha(1,6)-\betaGlcNAc (\phi,\psi)';
			T2 = 'Disaccharide: \alphaFucose-\alpha(1,6)-\betaGlcNAc (\omega,\psi)';
		case '16AB_MM'			% 2
			T  = 'Disaccharide: \alphaMannose-\alpha(1,6)-\betaMannose';
			T1 = 'Disaccharide: \alphaMannose-\alpha(1,6)-\betaMannose (\phi,\psi)';
			T2 = 'Disaccharide: \alphaMannose-\alpha(1,6)-\betaMannose (\omega,\psi)';
		case '16AB_MM_A'		% 3
			T  = 'Disaccharide: \alphaMannose-\alpha(1,6)-\alphaMannose';
			T1 = 'Disaccharide: \alphaMannose-\alpha(1,6)-\alphaMannose (\phi,\psi)';
			T2 = 'Disaccharide: \alphaMannose-\alpha(1,6)-\alphaMannose (\omega,\psi)';
		case '16BB_NM'			% 4
			T  = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaMannose';
			T1 = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaMannose (\phi,\psi)';
			T2 = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaMannose (\omega,\psi)';
		case '16BB_NG'			% 5
			T  = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaGalactose';
			T1 = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaGalactose (\phi,\psi)';
			T2 = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaGalactose (\omega,\psi)';
		case '16BB_NY'			% 6
			T  = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaGalNAc';
			T1 = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaGalNAc (\phi,\psi)';
			T2 = 'Disaccharide: \betaGlcNAc-\beta(1,6)-\betaGalNAc (\omega,\psi)';
	end %switch

	%% Initialization
	w		= M(1:L,4);
	wmax	= 360;    wmin  =    0;
	d1		= M(1:L,W - 3);	d2 = M(1:L,W - 2);
	fn1		= [dir,name,'.',num2str(seed),'.phi.eps'];		% Linkage plots: phi\psi
	fn2		= [dir,name,'.',num2str(seed),'.omega.eps'];	% Linkage plots: omega\psi
	%% Visualization
	figure(1)                          
	h12  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
	nx  = hist(x,xbin);
	bar(xbin,nx/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
	xlim([xmin xmax]);	ylabel('Frequency','FontSize',8);
	set(gca,'fontsize',7);	hold on;

	h13  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
	ny  = hist(y,ybin);
	bar(ybin,ny/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
	xlim([ymin ymax]);	view(h13,270,90);	set(gca,'xticklabel',[],'fontsize',7);
	ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);
	hold on;

	h14  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
	fastscatter(x,y,d1);	colorbar;	xlim([xmin xmax]);	ylim([ymin ymax]);	hold on;
	for i = 1 : size(dot,1)
        plot([dot(i,2) dot(i,2)], [ymin ymax],'--','color',color(i,:)); hold on;
        plot([xmin xmax], [dot(i,3) dot(i,3)],'--','color',color(i,:)); hold on;
        plot(dot(i,2),dot(i,3),'+','color',color(i,:));
    end % for
	grid on;	box on;
	xlabel('\phi(degree)','Fontsize',12);	suptitle(T1);		% No more input arguments allowed
	hold off;	print(fn1,'-depsc2');
	figure(2)
	h22  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
	nw  = hist(w,wbin);
	bar(wbin,nw/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
	xlim([wmin wmax]);	ylabel('Frequency','FontSize',8);
	set(gca,'fontsize',7);	hold on;

	h23  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
	ny  = hist(y,ybin);
	bar(ybin,ny/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
	xlim([ymin ymax]);	view(h23,270,90);	set(gca,'xticklabel',[],'fontsize',7);
	ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);	hold on;

	h24  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
	fastscatter(w,y,d2);	colorbar;	xlim([wmin wmax]);	ylim([ymin ymax]);	hold on;
	for i = 1 : size(dot,1)
        plot([dot(i,4) dot(i,4)], [ymin ymax],'--','color',color(i,:)); hold on;
        plot([wmin wmax], [dot(i,3) dot(i,3)],'--','color',color(i,:)); hold on;
        plot(dot(i,4),dot(i,3),'+','color',color(i,:));
    end % for
	grid on;	box on;	
	xlabel('\omega(degree)','Fontsize',12);	suptitle(T2); % No other input arguments allowed
	hold off;	print(fn2,'-depsc2');
	figure(3)
	fscatter3(x,y,w,d);		title(T);
	xlim([0 360]);	ylim([0 360]);	zlim([0 360]);
	xlabel('\phi(degree)','Fontsize',12);
	ylabel('\psi(degree)','Fontsize',12);
	zlabel('\omega(degree)','Fontsize',12);
	print(rn,'-depsc2');
end %if
end %function
