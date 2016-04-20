function [ oval ] = ovalIMG(file,name,ub,lb,stepsize)
% Use density info to get apporiate shape fit 
%
% INPUT VARIABLES:
%	file     - the density file
%                  1st  column - time
%                  2nd ~    column - phi,psi,omega
%                  last 2nd column - normalized density (0,1) 
%                  last     column - free energy
%                  * name without its directory
%   name    - the name of specific disaccharide
%       ub,lb    - upper bound/lower bound  
%	stepsize - density level increment
%       flag     - file format options(yet to be added):
%                  1 - plot w/o omega
%                  0 - plot w/  omega
% Copy-left: Cindy Lee,2014

%% Data loading and pre-processing
clc
dir    = strcat('/personae/cindy/Data/dimer/',name,'/');
if ischar(file)
    M  = load(strcat(dir,file));
elseif ismatrix(file)
	M  = file;
end % if
%% Initialization
J      = jet(50);
W      = size(M,2);
cut    = lb:stepsize:ub;
N      = length(cut) - 1;
oval   = zeros(length(cut) - 1, 8);

fh     = [dir,name,'.',num2str(lb),'~',num2str(ub),'.',num2str(N),'.eps'];
TT     = [strrep(name,'_','\_'),': Density limits = ',...
			num2str(lb),'~',num2str(ub),' with ',num2str(N),' ellipses'];
mh     = [dir,name,'.',num2str(lb),'~',num2str(ub),'.',num2str(N),'.dat'];
%% Visualization 
for i = 1 : N
	temp = M((M(:,W-1) > cut(i) & M(:,W-1) <= cut(i + 1)),:);		%increased density
	scatter(temp(:,2),temp(:,3),1,J(i,:),'o','filled'); daspect([1,1,1]);   hold on;
	ellipse = fit_ellipse(temp(:,2),temp(:,3));
	oval(i,1) = ellipse.X0_in;	oval(i,2) = ellipse.Y0_in;			%center for tilted ellipse
	oval(i,3) = ellipse.a;											%semi-major-axis
	oval(i,4) = ellipse.b;											%semi-minor-axis
	oval(i,5) = oval(i,3)/oval(i,4);								%ratio   
	oval(i,6) = ellipse.phi;                                        %ellipse tilt
	oval(i,7) = ellipse.X0;	oval(i,8) = ellipse.Y0;					%center for untilted ellipse
	set(findobj(gca,'type','line','color','k'),'color',J(i,:),'linewidth',1);	hold on;
end % for
hold off;	title(TT,'FontSize',14);	
xlabel('\phi','FontSize',14);	ylabel('\psi','FontSize',14);
print(fh,'-depsc2');
%% Data storage
%dlmwrite(mh,oval,'-append','delimiter','\t','newline','Unix','roffset',1);
dlmwrite(mh,oval,'delimiter','\t','newline','Unix');
end %function
