function [ ] = vmm_diff(name,vmm_s,D,K,N,varargin)
%% Energy difference between predicted model and simulation samples
%   Inputs:
%           name    -- Disaccharide name
%           vmm_s   -- Model predicted from VMM model 
%           D       -- Box size of simulation energy/pdf matrix
%           K       -- Number of components
%           N       -- Number of repeats for energy matrix
%           V1      -- C: Cutoff of energy/density for pdf comparison 
%           V2      -- S: Number of seeds in energy matrix
%   Copyright: Xindi Li (xindi.li@stonybrook.edu)
%% Initialization
R 	= -1.9872041 * 10^(-3);                            % Ideal gas constant
T 	= 300;                                                    % Temperature
C 	= 4.12;											  % >=100 dots in a box
% Parse inputs				
if nargin == 0
	return;
end
if nargin < 5
	error('TooFewInputs');
elseif nargin > 6
    S   = varargin{2};
    C   = varargin{1};
elseif nargin > 5
    C   = varargin{1};
end

% Results
dir     = strcat('~/Data/dimer/',name,'/');
dir1    = strcat(dir,'energy/');
Diff    = cell(N,360/D,360/D * K);               % Energy difference matrix
[XX,YY] = meshgrid(0:D:(360-D),0:D:(360-D));
Temp    = vmm_ang2rad([XX(:) YY(:)]);
%% Data Analysis
for i   = N
    e   = strcat(dir1,'name','.wt',num2str(D),'.',num2str(i),'_10.dat');
    Ene = flipud(load(e));
    Ene(Ene(:) == 100) = NaN;
    prob = exp(Ene/(R*T))/D^2;
    
    z   = pdf(Temp,vmm_s);
    z   = reshape(z,360/D,360/D);
    d   = prob - z/(sum(z(:))*D^2);
    Diff{r}(:,(k-1)*360/D + 1 : (k)*360/D) = d;
    
    dfile = strcat(dir,name,'_','Diff_',num2str(N),'_wt',num2str(D),'.txt');
    dlmwrite(dfile,Diff{i},'-append','delimiter','\t');
    
end
end % Function
