function [  ] = vmmtri_par(name,N,K1,K2,M,P1,P2,varargin)
%% VMMTRI_PAR fit trisaccharide glycosidic linkages using VMM model
% Initial guesses are based on parameters obtained in disaccharide cases
% and they will be given 5 degrees flexiblity
%   
%   Inputs:
%           Name -- Trisaccharide name consistent with the whole system
%           N    -- Repeated fit models with different initial guesses
%           K    -- Scalar : Largest component number allowed; or
%                   Vector : Chosen component numbers in model fittings
%           M    -- Sample ID from MC simulations
%           I    -- Initial guesses (could be a read-in file)
%           V1   -- s: Random number seed to generate initial guesses
%           V2   -- C: Cutoff for pdf comparison
%           V3   -- D: Resolution of simulation energy/pdf matrix
%
%   Note: This script is modified to enable parallel computing on a
%         multi-core cpu with MATLAB 2013b
%
%   Copyright: Xindi Li(xindi.li@stonybrook.edu)

%% Initialization
R 	= -1.9872041 * 10^(-3);                            % Ideal gas constant
T 	= 300;                                                    % Temperature
C 	= 4.12;											  % >=100 dots in a box
s   = 32; 						   % Random number seed for initial guesses
L   = 1;                                                % Resampling factor
% Parse inputs				
if nargin == 0
	return;
end
if nargin < 7
	error('TooFewInputs');
elseif nargin > 9
    s = varargin{1};
    C = varargin{2};
    D = varargin{3};
elseif nargin > 8
    s = varargin{1};
    C = vargrgin{2};
elseif nargin > 7
    s = varargin{1};
end

% [XX,YY] = meshgrid(0:D:(360-D),0:D:(360-D));
% Temp = vmm_ang2rad([XX(:) YY(:)]);

% Inital parameters
I1 = load(P1);
I2 = load(P2);
% Means                                 			% Randians in [-pi, pi]
rng(s);
Mu1 = cell(max(K1),N,2*max(K1));
Mu2 = cell(max(K2),N,2*max(K2));
for i = K1
    temp = I1(:,(i-1)*6 + 1: (i -1) *6 + 2);
    temp(~any(temp,2),:) = [];
    for j = 1 : i
        Mu1{i}(:, (j-1)*2+1 : 2*j) = [randi([floor(temp(j,1)-3) ...
            ceil(temp(j,1)+3)],N,1) ...
            randi([floor(temp(j,2)-3) ceil(temp(j,2)+3)],N,1)];
    end
end

for i = K2
    temp = I2(:,(i-1)*6+1: (i-1)*6+2);
    temp(~any(temp,2),:) = [];
    for j = 1 : i
        Mu2{i}(:, (j-1)*2+1 : 2*j) = [randi([floor(temp(j,1)-3) ...
            ceil(temp(j,1)+3)],N,1) ...
            randi([floor(temp(j,2)-3) ceil(temp(j,2)+3)],N,1)];
    end
end
% Kappa & Lambda
%kappa1 = cell(max(K1),
for i = K1
    kappa1(i)   = {I1(1:i,(i-1)*6 + 4 : (i-1)*6 + 5)};
    lambda1(i)  = {I1(1:i,i*6)};
end
for i = K2
    kappa2(i)   = {I2(1:i,(i-1)*6 + 4 : (i-1)*6 + 5)};
    lambda2(i)  = {I2(1:i,i*6)};
end

% Results
dir   = strcat('~/Data/trimer/',name,'/');
Table = cell(N,5,2*length(K1)*length(K2));            % Algorithm Performance
Param = cell(N,max(max(K1),max(K2))+1, 12*length(K1)*length(K2)); 
%% Data load
seed = load('~/Data/trimer/seed2.dat');

for j = M
    
    %f  = strcat(dir,name,'.1r30.1.dmatrix.5.dat');
    f   = strcat(dir,name,'.1.dmatrix.5.dat');
    TM = load(f);
    % Select every Lth point in radians [-pi, pi]
    A	= vmm_ang2rad(TM(1:L:end,2:3));
    B   = vmm_ang2rad(TM(1:L:end,4:5));
    clear TM;
    options = statset('MaxIter',500,'TolFun',1e-5);
    
    % matlabpool(3)
    for r = 1 : N
        Param{r} = zeros(max(max(K1),max(K2))+1,12*length(K1)*length(K2));
        % Loop for Non-reducing end angle fitting
        for k = K1
            id1 = find(K1 == k);
            S_nr = struct('Mu',vmm_ang2rad(reshape(reshape(Mu1{k}(r,:),...
                [1 2 k]),[2 k])'),'Kappa',kappa1{k},'Lambda',lambda1{k});
            tic
            vmm_nr = fitvmmdist(A,k,'Options',options,'Start',S_nr,...
              'Cortype','Sine');
            toc
        
            if vmm_nr.Converged
                
                for ii = 1 : length(K2)
                    Table{r}(:,(id1 -1)*2*length(K2)+(ii-1)*2+1) = ...
                    [vmm_nr.Ncomponents vmm_nr.AIC vmm_nr.BIC ...
                    vmm_nr.Converged vmm_nr.Iters]';
                    
                    Param{r}(1,(id1-1)*12*length(K2)+(ii-1)*12+1 : ...
                                (id1 -1)*12*length(K2) + (ii-1)*12+6) = ...
                                vmm_nr.Ncomponents;            
                    Param{r}(2:size(vmm_nr.Mu,1)+1, ...
                        (id1 -1)*12*length(K2) + (ii-1)*12+1 : ...
                        (id1 -1)*12*length(K2) + (ii-1)*12+6) = ...
                        [vmm_rad2ang(vmm_nr.Mu) vmm_nr.Pcomponents' ...
                        vmm_nr.Kappa vmm_nr.Lambda];                            
                end
                
                % Loop for reducing end angle fitting
                for l = K2
                    S_r = struct('Mu',vmm_ang2rad(reshape(reshape(Mu2{l}(r,:),...
                        [1,2,l]),[2 l])'),'Kappa',kappa2{l},'Lambda',lambda2{l});
                    tic
                    vmm_r = fitvmmdist(B,l,'Options',options,'Start',S_r,...
                        'Cortype','Sine');
                    toc
                    % Data collection
                     if vmm_r.Converged
                        id2 = find(K2 == l);
                        Table{r}(:,(id1-1)*2*length(K2) + 2*id2) = ...
                            [vmm_r.Ncomponents vmm_r.AIC vmm_r.BIC ...
                            vmm_r.Converged vmm_r.Iters]';
                        
                        Param{r}(1,(id1-1)*12*length(K2)+(id2-1)*12+7 : ...
                                (id1 -1)*12*length(K2) + id2*12) = ...
                                vmm_r.Ncomponents;  
                        Param{r}(2:size(vmm_r.Mu,1)+1, ...
                            (id1-1)*12*length(K2)+id2*12-5 : ...
                            (id1-1)*12*length(K2)+id2*12) = ...
                            [vmm_rad2ang(vmm_r.Mu) vmm_r.Pcomponents' ...
                            vmm_r.Kappa vmm_r.Lambda];        
                     else
                        Table{r}(:,(id1-1)*2*length(K2) + 2*id2) = ...
                        zeros(5,1);
                    
                        Param{r}(:,(id1-1)*12*length(K2)+id2*12-5 :...
                            (id1-1)*12*length(K2)+id2*12) = ...
                            zeros(max(max(K1),max(K2))+1,6);
                     end % Convergency 
                end % Reducing End Loop
            else
                for ii = 1 : length(K2)
                Table{r}(:,(id1 -1)*2*length(K2)+(ii-1)*2+1) = zeros(5,1);
                end 
            end
        end % Components
    end % Repeats
    %matlabpool close;
    
    % Write to files
    pfile = strcat(dir,name,'_','P_',num2str(N),'_',...
            num2str(seed(j)),'.txt');
    tfile = strcat(dir,name,'_','T_',num2str(N),'_',...
                num2str(seed(j)),'.txt');
    %dfile = strcat(dir,name,'_','D_',num2str(N),'R_',...
%                num2str(seed(j)),'.txt');
            
    for ii = 1 :N
        dlmwrite(tfile,Table{ii},'-append','delimiter','\t');
        dlmwrite(pfile,Param{ii},'-append','delimiter','\t');
        %dlmwrite(dfile,Diff{ii},'-append','delimiter','\t');
    end
end % Seeds
end % Function

