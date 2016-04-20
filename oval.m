function [ oval ] = oval(file,name,seed,ub,lb,stepsize)
% Use density info to get apporiate shape fit 
%
% INPUT VARIABLES:
%	file     - the density file
%                  1st  column - time
%                  2nd ~    column - phi,psi,omega
%                  last 2nd column - normalized density (0,1) 
%                  last     column - free energy
%                  * name without its directory
%   name     - the name of specific disaccharide
%   seed     - generating seed of one sample
%   ub,lb    - upper bound/lower bound  
%	stepsize - density level increment
%       flag     - file format options(yet to be added):
%                  1 - plot w/o omega
%                  0 - plot w/  omega
% NOTE : For display of found ellipse, please use ovalIMG 
% NOTE : Data would be accumulated 
% Copy-left: Cindy Lee,2014

%% Data loading and pre-processing
clc
dir    = strcat('/personae/cindy/Data/dimer/',name,'/density/');
if ischar(file)
    M  = load(strcat(dir,file));
elseif ismatrix(file)
	M  = file;
end % if
%% Initialization
W      = size(M,2);
cut    = lb:stepsize:ub;
N      = length(cut) - 1;
oval   = zeros(length(cut) - 1, 8);
mh     = [dir,name,'.',num2str(seed),'.',num2str(lb),'~',num2str(ub),'.',num2str(N),'.dat'];
%% Visualization 
for i = 1 : N
    if  cut(i + 1) <= 0.98
        temp = M((M(:,W-1) > cut(i) & M(:,W-1) <= cut(i + 1)),:);		%increased density
        ellipse = fit_ellipse(temp(:,2),temp(:,3));
        oval(i,1) = ellipse.X0_in;	oval(i,2) = ellipse.Y0_in;			%center for tilted ellipse
        oval(i,3) = ellipse.a;											%semi-major-axis
        oval(i,4) = ellipse.b;											%semi-minor-axis
        oval(i,5) = oval(i,3)/oval(i,4);								%ratio   
        oval(i,6) = ellipse.phi;                                        %ellipse tilt
        oval(i,7) = ellipse.X0;	oval(i,8) = ellipse.Y0;					%center for untilted ellipse
    elseif cut(i + 1) > 0.98 && cut(i + 1) ~= 1 
        temp = M(M(:,W-1) > cut(i),:);
        ellipse = fit_ellipse(temp(:,2),temp(:,3));
        oval(i,1) = ellipse.X0_in;	oval(i,2) = ellipse.Y0_in;			%center for tilted ellipse
        oval(i,3) = ellipse.a;											%semi-major-axis
        oval(i,4) = ellipse.b;											%semi-minor-axis
        oval(i,5) = oval(i,3)/oval(i,4);								%ratio   
        oval(i,6) = ellipse.phi;                                        %ellipse tilt
        oval(i,7) = ellipse.X0;	oval(i,8) = ellipse.Y0;					%center for untilted ellipse
    else
        temp = M(M(:,W-1) > cut(i-1),:);
        ellipse = fit_ellipse(temp(:,2),temp(:,3));
        oval(i,1) = ellipse.X0_in;	oval(i,2) = ellipse.Y0_in;			%center for tilted ellipse
        oval(i,3) = ellipse.a;											%semi-major-axis
        oval(i,4) = ellipse.b;											%semi-minor-axis
        oval(i,5) = oval(i,3)/oval(i,4);								%ratio   
        oval(i,6) = ellipse.phi;                                        %ellipse tilt
        oval(i,7) = ellipse.X0;	oval(i,8) = ellipse.Y0;					%center for untilted ellipse
    end
end % for
%% Data storage
%dlmwrite(mh,oval,'-append','delimiter','\t','newline','Unix','roffset',1);
dlmwrite(mh,oval,'delimiter','\t','newline','Unix');
end % function
