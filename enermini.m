function [ ] = enermini(file,S,stepsize,sim,flag)
% energy matrix  --> Positions of local minima
% 
% INPUT VARIABLES:
%       file     - the rescaled linkage data whose
%                  1st	column - time
%                  2nd ~    column - phi,psi,omega
%                  last     column - free energy
%                  * name without its directory
%       S        - the number of seeds used to generate energy landscapes 
%                  ( * 10^4)
%       stepsize - the length of the energy box
%       sim      - string for specific seed combination 
%		flag	 - 1: linkage w/o extra omega
%				 - 0: linkage w/  extra omega 
% OUTPUT:
%       #column  1st               2nd          3rd           4th  5th  6th
%                linear index      line index   column index  x    y    energy
% Copy-left: Cindy Lee,2014

clc
str		= regexp(file,'\.','split');
name	= char(str(1));
dir		= strcat('~/Data/dimer/',name,'/energy/');
M		= load(strcat(dir,file));
dim		= size(M,2);
xmax	= 360;	xmin = 0;
ymax	= 360;	ymin = 0;
zmax	= 360;	zmin = 0;
mini	= [dir,num2str(S),'s.',num2str(stepsize),'.',sim,'.mini.mat'];
%% Computation
if flag 
	[~,~,XMIN,IMIN] = extrema2(M);
	m  = size(IMIN,1);
	I  = zeros(m,2);      % to store the according (i,j) position
	F  = zeros(m,2);      % to store the according energy box 
	for i = 1 : m
    	if  mod(IMIN(i),size(M,1)) == 0
            I(i,1) = size(M,1);
	    	I(i,2) = IMIN(i)/size(M,1);
		else
			I(i,2) = floor(IMIN(i)/size(M,1)) + 1;
			I(i,1) = mod(IMIN(i),size(M,1)); 
    	end %if
	end %for
	% Get real positions
	for j = 1 : m
    	F(j,1) = (I(j,2) - 1) * stepsize + xmin;			%phi 
    	F(j,2) = ymax  - I(j,1) * stepsize;					%psi
	end %for
	% Data Collection  
	switch S
		case 100 
		data100 = [IMIN I F XMIN]; save(mini,'data100');
		case 70 
		data70 = [IMIN I F XMIN]; save(mini,'data70');
		case 50
		data50 = [IMIN I F XMIN]; save(mini,'data50');
		case 30
		data30 = [IMIN I F XMIN]; save(mini,'data30');
		case 20
		data20 = [IMIN I F XMIN]; save(mini,'data20');
		case 10
		data10 = [IMIN I F XMIN]; save(mini,'data10');
	end % switch 
else
	% Convert 2D energy matrix to according 3D format
	M3 = zeros(dim,dim,dim);
	for k = 1 : dim
		M3(:,:,k) =  M( ((k - 1) * dim + 1 ) : k * dim, :);	
	end
	[~,~,XMIN,IMIN] = MinimaMaxima3D(M3,1,1);
	m  = size(IMIN,1);
	F  = zeros(m,3);      % to store the real energy box 
	% Get real positions
	for j = 1 : m
    	F(j,1) = (IMIN(j,2) - 1) * stepsize + xmin;				%phi 
    	F(j,2) = ymax  - IMIN(j,1) * stepsize;					%psi
		F(j,3) = (IMIN(j,3) - 1) * stepsize + zmin;
	end %for
    % Data Collection  
    switch S
        case 100
        data100 = [IMIN F XMIN]; save(mini,'data100');
        case 70
        data70 = [IMIN F XMIN]; save(mini,'data70');
        case 50
        data50 = [IMIN F XMIN]; save(mini,'data50');
        case 30
        data30 = [IMIN F XMIN]; save(mini,'data30');
        case 20
        data20 = [IMIN F XMIN]; save(mini,'data20');
        case 10
        data10 = [IMIN F XMIN]; save(mini,'data10');
    end % switch 
end % function 
