function [ ] = glylinkage_tri(file,S,type)
% Glycosidic linkage --> Ramachandran plot (phi,psi,omega) for DISACCHARIDES
% Scatter plot with color indicating data density
% INPUT VARIABLES:
%       file  - is the density matrix file
%               1st	column - time
%               2nd ~    column - phi,psi,omega
%               2nd last column - normalized density
%               last     column - free energy
%               * name without its directory
%       S     - is the number of seeds used to generate the plot
%       type  - file format options:
%               11 - plot w/o omega
%               10 - only 2nd plot w/ omega
%				01 - only 1st plot w/ omega
%				00 - both w/ omega
% Copy-left: Cindy Lee,2014
clc
str		= regexp(file,'\.','split'); name = char(str(1));		% Trisaccharide Name
dir		= strcat('~/Data/trimer/',name,'/');					% File Directory
lh		= [dir,name,'.',num2str(S),'.lhist.dat'];				% Linkage Linear Histogram
dh		= [dir,name,'.',num2str(S),'.dhist.dat'];				% Density Histogram
M		= load(strcat(dir,file));
W		= size(M,2);	L = S * 1e05;							% Data size
dot		= M((M(:,W) == min(M(:,W))),:);							% Lowest energy matrix
color	= hsv(10);	dbin	= 0 : 0.05 : 0.95;					% Density bins
xmax	= 360;	xmin = 0;	xbin = (xmin + 5) :10: (xmax - 5);	% Linkage bins	
ymax	= 360;	ymin = 0;	ybin = (ymin + 5) :10: (ymax - 5);
zmax	= 360;	zmin = 0;	zbin = (zmin + 5) :10: (zmax - 5);

fn1		= [dir,name,'.',num2str(S),'seeds.n.eps'];					% for linkage plots: non-reducing end
fr2		= [dir,name,'.',num2str(S),'seeds.r.eps'];					% for linkage plots: reducing end

% Load specific data
switch type
	case '11' 
		x	= M(1:L,2);	y	= M(1:L,3);							% phi/psi: link 1
		u	= M(1:L,4);	v	= M(1:L,5);							% phi/psi: link 2
		d1	= M(1:L,6); d2	= M(1:L,7);							% Density
	case '10'
		x   = M(1:L,2); y   = M(1:L,3);							% phi/psi    : link 1
		u   = M(1:L,4); v   = M(1:L,5);	w	= M(1:L,6);			% phi/psi/chi: link 2
		d1  = M(1:L,7);	d21	= M(1:L,8);	d22	= M(1:L,9);         % Density
	case '01'
		x   = M(1:L,2); y   = M(1:L,3); z	= M(1:L,4);			% phi/psi/chi: link 1
		u   = M(1:L,5); v   = M(1:L,6);							% phi/psi    : link 2
		d11 = M(1:L,7); d12	= M(1:L,8);	d2	= M(1:L,9);
	case '00'
		x   = M(1:L,2); y   = M(1:L,3); z   = M(1:L,4);         % phi/psi/chi: link 1
		u   = M(1:L,5); v   = M(1:L,6);	w	= M(1:L,7);			% phi/psi/chi: link 2
		d11	= M(1:L,8);	d12	= M(1:L,9);
		d21 = M(1:L,10);d22 = M(1:L,11);
end % Switch type
%% Specific Data Sets Parameters
switch type
	case '11' 
		switch name
			case '22MMM'
                T = 'Trisaccharide: \alphaMannose-(1,2)-\alphaMannose-(1,2)-\alphaMannose';
                pos1x= 100;	pos1y = 100;	pos2x= 100;	pos2y = 100;
                label1 = 'Non-Reducing end : \alpha(1,2) linkage' ;
                label2 = 'Reducing end: \alpha(1,2) linkage ';
			case '23MMM'
                T = 'Trisaccharide: \alphaMannose-(1,2)-\alphaMannose-(1,3)-\alphaMannose';
                pos1x= 200;	pos1y = 150;	pos2x= 200;	pos2y = 150;
                label1 = 'Non-reducing end: \alpha(1,2) linkage';
                label2 = 'Reducing end: \alpha(1,3) linkage';
			case '42GNM'
                T = 'Trisaccharide: \betaGalactose-(1,4)-\betaGlcNAc-(1,2)-\alphaMannose';
                pos1x= 200; pos1y = 150;	pos2x= 200; pos2y = 150;
                label1 = 'Non-reducing end: \beta(1,4) linkage';
                label2 = 'Reducing end : \beta(1,2) linkage' ;
			case '43GNG'
                T = ' Trisaccharide: \betaGalactose-(1,4)-\betaGlcNAc-(1,3)-\betaGalactose';
                pos1x= 200; pos1y = 150;	pos2x= 200; pos2y = 150;
                label1 = 'Non-reducing end: \beta(1,4) linkage';
                label2 = 'Reducing end : \beta(1,3) linkage' ;
			case '44MNN'
                T = 'Trisaccharide: \betaMannose-(1,4)-\betaGlcNAc-(1,4)-\betaGlcNAc';
                pos1x = 200; pos1y = 150;	pos2x= 200; pos2y = 150;
                label1 = 'Non-reducing end: \beta(1,4) linkage';
                label2 = 'Reducing end : \beta(1,4) linkage' ;
            case '24NMN'
                T = 'Trisaccharide: \betaGlcNAc -(1,2) - \alphaMannose - (4,1)-\betaGlcNAc';
                pos1x = 200; pos1y = 150;	pos2x= 200; pos2y = 150;
                label1 = '\beta(1,2) linkage';
                label2 = '\beta(1,4) linkage';
            case '44GNM'
                T = 'Trisaccharide: \betaGalactose-(1,4)-\betaGlcNAc-(1,4)-\betaMannose';
                pos1x = 75; pos1y = 200;	pos2x= 75; pos2y = 200;
                label1 = 'Non-reducing end: \beta(1,4) linkage';
                label2 = 'Reducing end : \beta(1,4) linkage' ;
            case '43YNF'
                T = 'Trisaccharide:\betaGalNAc-(1,4)-\betaGlcNAc-(3,1)-\alphafucose';
                pos1x = 25; pos1y = 250;    pos2x = 75; pos2y = 100;
                label1 = 'Non-reducing end : \beta(1,4) linkage';
                label2 = 'Fucose attachment : \alpha(1,3) linkage';
			case '43YGS'
				T = 'Trisaccharide:\betaGalNAc-(1,4)-\betaGalactose-(3,2)-\alphaNeu5Ac';
				pos1x = 25; pos1y = 250;	pos2x = 75; pos2y = 100;
				label1 = 'Non-reducing end :\beta(1,4) linkage';
				label2 = 'Fucose attachment : \alpha(3,2) linkage';
			case '32YGF'
				T = 'Trisaccharide:\alphaGalNAc-(1,3)-\betaGalactose-(2,1)-\alphaFucose';
				pos1x = 25; pos1y = 250;	pos2x = 75; pos2y = 100;
				label1 = 'Non-reducing end :\alpha(1,3) linkage';
				label2 = 'Fucose attachment : \alpha(1,2) linkage';
			case '32GGF'
				T = 'Trisaccharide:\alphaGalactose-(1,3)-\betaGalactose-(2,1)-\alphaFucose';
				pos1x = 25; pos1y = 250;	pos2x = 75; pos2y = 100;
				label1 = 'Non-reducing end :\alpha(1,3) linkage';
				label2 = 'Fucose attachment : \alpha(1,2) linkage';
		end % Switch Name
	case '10'
		switch name
			case '26MMM'
                T   = 'Trisaccharide: \alphaMannose-(1,2)-\alphaMannose-(1,6)-\alphaMannose';
                pos1x = 100;	pos1y = 100;	
                pos2x = 200;	pos2y = 300;	pos2z = 200;
                label1 = 'Non-reducing end: \alpha(1,2) linkage';
                label2 = 'Reducing end: \alpha(1,2) linkage';
			case '36MMM_B'
                T = 'Trisaccharide: \alphaMannose-(1,2)-\betaMannose-(6,1)-\alphaMannose';
                pos1x = 200;	pos1y = 100;
                pos2x = 200;	pos2y = 300;	pos2z = 200;
                label1 = '\alpha(1,3) linkage';
                label2 = '\alpha(1,6) linkage';
			case '36MMM_C'
                T = 'Trisaccharide: \alphaMannose-(1,2)-\betaMannose-(6,1)-\alphaMannose';
                pos1x = 200;	pos1y = 100;
                pos2x = 200;	pos2y = 300;	pos2z = 200;
                label1 = '\alpha(1,3) linkage';
                label2 = '\alpha(1,6) linkage';
            case '46GNM'
                T = 'Trisaccharide: \betaGal-(1,4)-\betaGlcNAc-(1,6)-\alphaMannose';
                pos1x = 75;	pos1y = 200;
                pos2x = 200;	pos2y = 300;	pos2z = 200;
                label1 = 'Non-reducing end: \alpha(1,4) linkage';
                label2 = 'Reducing end: \beta(1,6) linkage';
            case '26NMN'
                T = 'Trisaccharide: \betaGlcNAc-(1,2)-\alphaMannose-(6,1)-\betaGlcNAc';
                pos1x = 75;	pos1y = 150;
                pos2x = 50;	pos2y = 325;	pos2z = 150;
                label1 = '\beta(1,2) linkage';
                label2 = '\beta(1,6) linkage';
            case '36NGN_C'
                T = 'Trisaccharide: \betaGlcNAc-(1,3)-\betaGal-(6,1)-\betaGlcNAc';
                pos1x = 75;	pos1y = 100;
                pos2x = 50;	pos2y = 325;	pos2z = 150;
                label1 = '\beta(1,3) linkage';
                label2 = '\beta(1,6) linkage';
            case '36GYN_C'
                T = 'Trisaccharide: \betaGal-(1,3)-\alphaGalNAc-(6,1)-\betaGlcNAc';
                pos1x = 75;	pos1y = 100;
                pos2x = 50;	pos2y = 325;	pos2z = 150;
                label1 = '\beta(1,3) linkage';
                label2 = '\beta(1,6) linkage';
            case '46NNF'
                T = 'Trisaccharide: \betaGlcNAc-(1,4)-\betaGlcNAc-(6,1)-\alphaFucose';
                pos1x = 75;	pos1y = 100;
                pos2x = 50;	pos2y = 325;	pos2z = 150;
                label1 = '\beta(1,4) linkage';
                label2 = '\alpha(1,6) linkage';
		end % Switch Name
	case '01'
	case '00'
end % Switch Type
%% Visualization
%Ramachandron Plot of non-reducing end linkage
switch type
	case {'11','10'}
		figure(1) 
		h12  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
		nx  = hist(x,xbin);
		bar(xbin,nx/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([xmin xmax]);
		ylabel('Frequency','FontSize',8);	set(gca,'fontsize',7);	hold on;

		h13  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
		ny  = hist(y,ybin);
		bar(ybin,ny/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([ymin,ymax]);	view(h13,270,90);	set(gca,'xticklabel',[],'fontsize',7);
		ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);

		h14  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
		fastscatter(x,y,d1);	colorbar;	xlim([xmin xmax]); ylim([ymin ymax]);   hold on;
		for i = 1 : size(dot,1)
		plot([dot(i,2) dot(i,2)], [ymin ymax],'--','color',color(i,:));
		plot([xmin xmax], [dot(i,3) dot(i,3)],'--','color',color(i,:));
        plot(dot(i,2),dot(i,3),'+','color',color(i,:));
		end % for
		grid on;	box on;
		text(pos1x,pos1y,label1,'fontsize',8);
		xlabel('\phi(degree)','Fontsize',12);	suptitle(T);    % No other input arguments allowed
		hold off;	print(fn1,'-depsc2');
	case {'01','00'}
		fn11	= [dir,name,'.',num2str(S),'.n.phi.eps'];		% Linkage plots: phi\psi
		fn12	= [dir,name,'.',num2str(S),'.n.omega.eps'];	% Linkage plots: omega\psi
		
		figure(11)
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

		h14  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
		fastscatter(x,y,d11);	colorbar;	xlim([xmin xmax]);	ylim([ymin ymax]);  hold on;
		for i = 1 : size(dot,1)
        plot([dot(i,2) dot(i,2)], [ymin ymax],'--','color',color(i,:)); 
        plot([xmin xmax], [dot(i,3) dot(i,3)],'--','color',color(i,:)); 
        plot(dot(i,2),dot(i,3),'+','color',color(i,:));
        end % for
		grid on;	box on;	text(pos1x,pos1y,label1,'fontsize',8);
		xlabel('\phi(degree)','Fontsize',12);	suptitle(T);		% No more input arguments allowed
		hold off;	print(fn11,'-depsc2');

		figure(12)
		h22  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
		nz  = hist(z,zbin);
		bar(zbin,nz/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([zmin zmax]);	ylabel('Frequency','FontSize',8);	set(gca,'fontsize',7);	hold on;
		set(gca,'fontsize',7);	hold on;

		h23  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
		ny  = hist(y,ybin);
		bar(ybin,ny/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([ymin ymax]);	view(h23,270,90);	set(gca,'xticklabel',[],'fontsize',7);
		ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);

		h24  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
		fastscatter(z,y,d12);	colorbar;	xlim([zmin zmax]);	ylim([ymin ymax]);
		for i = 1 : size(dot,1)
        	plot([dot(i,4) dot(i,4)], [ymin ymax],'--','color',color(i,:)); hold on;
        	plot([zmin zmax], [dot(i,3) dot(i,3)],'--','color',color(i,:)); hold on;
        	plot(dot(i,4),dot(i,3),'+','color',color(i,:));
    	end % for
		grid on;	box on;	text(pos1z,pos1y,label1,'fontsize',8);
		xlabel('\omega(degree)','Fontsize',12);	suptitle(T);	% No other input arguments allowed
		hold off;	print(fn1,'-depsc2');
end % switch
%Ramachandron Plot of reducing end linkage
switch type
	case {'11','01'}
		figure(2) 
		h32  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
		nu  = hist(u,xbin);
		bar(xbin,nu/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([xmin xmax]);	ylabel('Frequency','FontSize',8);	set(gca,'fontsize',7);	hold on;

		h33  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
		nv  = hist(v,ybin);
		bar(ybin,nv/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([ymin,ymax]);	view(h33,270,90);	set(gca,'xticklabel',[],'fontsize',7);
		ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);

		h34  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
		fastscatter(u,v,d2);	colorbar;	xlim([xmin xmax]); ylim([ymin ymax]);	hold on;	
		if size(dot,2) == 8
			for i = 1 : size(dot,1)
			plot([dot(i,4) dot(i,4)], [ymin ymax],'--','color',color(i,:));
			plot([xmin xmax], [dot(i,5) dot(i,5)],'--','color',color(i,:));
        	plot(dot(i,4),dot(i,5),'+','color',color(i,:));
			end % for
		else
			for i = 1 : size(dot,1)
            plot([dot(i,5) dot(i,5)], [ymin ymax],'--','color',color(i,:));
            plot([xmin xmax], [dot(i,6) dot(i,6)],'--','color',color(i,:));
            plot(dot(i,5),dot(i,6),'+','color',color(i,:));
            end % for
		end % if 
		grid on;	box on;	text(pos2x,pos2y,label2,'fontsize',8);
		xlabel('\phi(degree)','Fontsize',12);	suptitle(T);    % No other input arguments allowed
		hold off;	print(fr2,'-depsc2');
	case {'10','00'}
		fr21	= [dir,name,'.',num2str(S),'.r.phi.eps'];		% Linkage plots: phi\psi
		fr22	= [dir,name,'.',num2str(S),'.r.omega.eps'];		% Linkage plots: omega\psi
		
		figure(21)
		h42  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
		nu  = hist(u,xbin);
		bar(xbin,nu/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([xmin xmax]);	ylabel('Frequency','FontSize',8);
		set(gca,'fontsize',7);	hold on;

		h43  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
		nv  = hist(v,ybin);
		bar(ybin,nv/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([ymin ymax]);	view(h43,270,90);	set(gca,'xticklabel',[],'fontsize',7);
		ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);

		h44  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
		fastscatter(u,v,d21);	colorbar;	xlim([xmin xmax]);	ylim([ymin ymax]);  hold on;
        if size(dot,2) == 12
			for i = 1 : size(dot,1)
			plot([dot(i,5) dot(i,5)], [ymin ymax],'--','color',color(i,:));
			plot([xmin xmax], [dot(i,6) dot(i,6)],'--','color',color(i,:));
        	plot(dot(i,5),dot(i,6),'+','color',color(i,:));
			end % for
		else
			for i = 1 : size(dot,1)
            plot([dot(i,4) dot(i,4)], [ymin ymax],'--','color',color(i,:));
            plot([xmin xmax], [dot(i,5) dot(i,5)],'--','color',color(i,:));
            plot(dot(i,4),dot(i,5),'+','color',color(i,:));
            end % for
		end % if 
		grid on;	box on;	text(pos2x,pos2y,label2,'fontsize',8);
		xlabel('\phi(degree)','Fontsize',12);	suptitle(T);		% No more input arguments allowed
		hold off;	print(fr21,'-depsc2');

		figure(22)
		h52  = subplot(2,2,2,'Position',[0.3,0.75,0.484,0.1]);
		nw  = hist(w,zbin);
		bar(zbin,nw/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([zmin zmax]);	ylabel('Frequency','FontSize',8);	set(gca,'fontsize',7);	hold on;
		set(gca,'fontsize',7);	hold on;

		h53  = subplot(2,2,3,'Position',[0.16,0.1,0.08,0.6]);
		nv  = hist(v,ybin);
		bar(ybin,nv/L,'FaceColor',[0 .5 .5],'EdgeColor',[0.5 0.5 0.5],'barwidth',0.9);
		xlim([ymin ymax]);	view(h53,270,90);	set(gca,'xticklabel',[],'fontsize',7);
		ylabel('Frequency','FontSize',8);	xlabel('\psi(degree)','FontSize',12);

		h54  = subplot(2,2,4,'Position',[0.3,0.1,0.6,0.6]);
		fastscatter(w,v,d22);	colorbar;	xlim([zmin zmax]);	ylim([ymin ymax]);  hold on;
        if size(dot,2) == 12
			for i = 1 : size(dot,1)
			plot([dot(i,7) dot(i,7)], [ymin ymax],'--','color',color(i,:));
			plot([xmin xmax], [dot(i,6) dot(i,6)],'--','color',color(i,:));
        	plot(dot(i,7),dot(i,6),'+','color',color(i,:));
			end % for
		else
			for i = 1 : size(dot,1)
            plot([dot(i,6) dot(i,6)], [ymin ymax],'--','color',color(i,:));
            plot([xmin xmax], [dot(i,5) dot(i,5)],'--','color',color(i,:));
            plot(dot(i,7),dot(i,5),'+','color',color(i,:));
            end % for
		end % if 
		grid on;	box on;	text(pos2z,pos2y,label2,'fontsize',8);
		xlabel('\omega(degree)','Fontsize',12);	suptitle(T); % No other input arguments allowed
		hold off;	print(fr22,'-depsc2');
end % switch
%% Data Collection
switch type
    case '11'
	denh1 = hist(d1,dbin);
	denh2 = hist(d2,dbin);
	data2 = [xbin' nx'/L ny'/L nu'/L nv'/L];
	data3 = [dbin' denh1'/L denh2'/L];
	dlmwrite(lh,data2,'delimiter','\t','newline','Unix');
	dlmwrite(dh,data3,'delimiter','\t','newline','Unix');
    case '10'
	denh1  = hist(d1,dbin);
    denh21 = hist(d21,dbin); denh22 = hist(d22,dbin);
    data2  = [xbin' nx'/L  ny'/L nu'/L nv'/L nw'/L];
    data3  = [dbin' denh1'/L denh21'/L denh22'/L];
    dlmwrite(lh,data2,'delimiter','\t','newline','Unix');
    dlmwrite(dh,data3,'delimiter','\t','newline','Unix');
    case '01'
	denh11 = hist(d11,dbin); denh12 = hist(d12,dbin);
    denh2  = hist(d21,dbin); 
    data2  = [xbin' nx'/L  ny'/L nz'/L nu'/L nv'/L ];
    data3  = [dbin' denh11'/L denh12'/L denh2'/L];
    dlmwrite(lh,data2,'delimiter','\t','newline','Unix');
    dlmwrite(dh,data3,'delimiter','\t','newline','Unix');
	case '00'
	denh11 = hist(d11,dbin); denh12 = hist(d12,dbin);
	denh21 = hist(d21,dbin); denh22 = hist(d22,dbin);
	data2  = [xbin' nx'/L  ny'/L nz'/L nu'/L nv'/L nw'/L];
	data3  = [dbin' denh11'/L denh12'/L denh21'/L denh22'/L];
    dlmwrite(lh,data2,'delimiter','\t','newline','Unix');
    dlmwrite(dh,data3,'delimiter','\t','newline','Unix');
end % Switch type
end %function
