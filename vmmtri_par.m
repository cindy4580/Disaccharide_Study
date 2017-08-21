function [  ] = vmmtri_par(name,N,K1,K2,M,P1,P2,type,varargin)
%% VMMTRI_PAR fitting starter script
% Initial guesses are based on parameters obtained in disaccharide cases
% and the angular mean are within 10 degrees away; component weights are
% randomly generated
%   
%   Inputs:
%           Name -- Trisaccharide name consistent with the whole system
%           N    -- Repeated fitting with different initial guesses
%           K    -- Number of components at most
%           M    -- Sample ID from MC simulations
%           P1   -- Non-reducing end initial guesses
%           P2   -- Reducing end initial guesses
%           Type -- File format options
%                   11 - two linkages w/o omega
%                   10 - the second linkage w/ omega
%           V1   -- s: Random number seed to generate initial guesses
%           V2   -- L: Resampling factor
%
%   Note: This script is modified to enable parallel computing on a
%         multi-core cpu with MATLAB 2013b therefore contains outdated
%         built-in functions
%
%   Copyright: Xindi Li(xindi.li@stonybrook.edu)

%% Initialization
L = 1;                                                  % Resampling factor
s = 32; 						   % Random number seed for initial guesses
% Parse inputs				
if nargin == 0
	return;
end
if nargin < 8
	error('TooFewInputs');
elseif nargin > 9
    s = varargin{1};
    L = varargin{2};
elseif nargin > 8
    s = varargin{1};
end
% Inital parameters
rng(s);
I1  = load(P1);
I2  = load(P2);
in1 = 6;
if strcmp(type,'11')
    in2 = 6;
elseif strcmp(type,'10')
    in2 = 10;
end
Mu1 = cell(max(K1),N,3*max(K1));
Mu2 = cell(max(K2),N,3*max(K2));
Ka1 = cell(max(K1),N,3*max(K1));
Ka2 = cell(max(K2),N,3*max(K2));
L1  = cell(max(K1),N,3*max(K1));
L2  = cell(max(K2),N,3*max(K2));
% Mean
for i = K1
    temp = I1(1:i ,(i-1)*in1+1 : (i-1)*in1+2);
    if any(temp)
        for j = 1 : i
            Mu1{i}(:, (j-1)*2+1 : 2*j) = [randi([floor(temp(j,1)-5) ...
                ceil(temp(j,1)+5)],N,1) ...
                randi([floor(temp(j,2)-5) ceil(temp(j,2)+5)],N,1)];
        end
    end
end

if strcmp(type,'11')
    for i = K2
        temp = I2(1:i, (i-1)*in1+1 : (i-1)*in1+2);
        if any(temp)
            for j = 1 : i
                Mu2{i}(:, (j-1)*2+1 : 2*j) = [randi([floor(temp(j,1)-5) ...
                    ceil(temp(j,1)+5)],N,1) ...
                    randi([floor(temp(j,2)-5) ceil(temp(j,2)+5)],N,1)];
            end
        end
    end
elseif strcmp(type,'10')
    for i = K2
        temp = I2(1:i, (i-1)*in2+1 : (i -1)*in2+3);
        if any(temp)
            for j = 1 : i
                Mu2{i}(:, (j-1)*3+1 : 3*j) = [randi([floor(temp(j,1)-5) ...
                    ceil(temp(j,1)+5)],N,1) ...
                    randi([floor(temp(j,2)-5) ceil(temp(j,2)+5)],N,1) ...
                    randi([floor(temp(j,3)-5) ceil(temp(j,3)+5)],N,1)];
            end
        end
    end
else
    warning('Unknown linkage type at the reducing end');
end
% Kappa & Lambda
for i = K1
    Ka1{i}  = I1(1:i,(i-1)*in1 + 4 : (i-1)*in1 + 5);
    L1{i}   = I1(1:i,i*in1);
end
if strcmp(type,'11')
    for i = K2
        Ka2{i} = I2(1:i,(i-1)*in1 + 4 : (i-1)*in1 + 5);
        L2{i}  = I2(1:i,i*in1);
    end
elseif strcmp(type,'10')
    for i = K2
        Ka2{i} = I2(1:i,(i-1)*in2 + 5 : (i-1)*in2 + 7);
        L2{i}  = I2(1:i,(i-1)*in2 + 8 : (i-1)*in2 + 10);
    end
else
    warning('Unknown linkage type at the reducing end');
end

% Results
dir   = strcat('~/Data/trimer/',name,'/');
Table = cell(N,5,2*length(K1)*length(K2));          % Algorithm Performance
Param = cell(N,max(max(K1),max(K2))+1, in1*max(K1)+in2*max(K2)); 
%% Data load
sd      = load('~/Data/dimer/seed1.dat');
for j   = M
    f   = strcat(dir,'density/',name,'.',num2str(sd(j)),'.d.dat');
    TM  = load(f);
    % Select every Lth point in radians [-pi, pi]
    A	= vmm_ang2rad(TM(1:L:end,2:3));   
    if strcmp(type,'11')
        B   = vmm_ang2rad(TM(1:L:end,4:5));
    elseif strcmp(type,'10')
        B   = vmm_ang2rad(TM(1:L:end,4:6));
    else
        warning('Unknown Linkage type at the reducing end');
    end    
    clear TM;
    options = statset('MaxIter',500,'TolFun',1e-5);
    
    % matlabpool(3)
    for r = 1: N
        Param{r} = zeros(max(max(K1),max(K2))+1,in1*max(K1)+in2*max(K2));
        % Loop for Non-reducing end angle fitting
%         for k = K1
%             if any(Mu1{k})
%                 S_nr = struct('Mu',vmm_ang2rad(reshape(reshape(Mu1{k}(r,:),...
%                     [1 2 k]),[2 k])'),'Kappa',Ka1{k},'Lambda',L1{k});
%                 tic
%                 vmm_nr = fitvmmdist(A,k,'Options',options,'Start',S_nr,...
%                   'Cortype','Sine');
%                 toc
%             else
%                 error('No convergent initial guesses in disaccharides');
%             end
%             Param{r}(1,in1*(k-1) + 1 : in1*k) = vmm_nr.Ncomponents;  
%             if vmm_nr.Converged
%                 Table{r}(:,k)= [vmm_nr.Ncomponents vmm_nr.AIC ...
%                                 vmm_nr.BIC vmm_nr.Converged vmm_nr.Iters]';
%                     
%                 Param{r}(2:size(vmm_nr.Mu,1)+1, (k-1)*in1+1 : k*in1) = ...         
%                         [vmm_rad2ang(vmm_nr.Mu) vmm_nr.Pcomponents' ...
%                          vmm_nr.Kappa vmm_nr.Lambda];                            
%             else % NR not converged
%                 Table{r}(:,k)= NaN(5,1);
%                 Param{r}(2:max(max(K1),max(K2))+1, (k-1)*in1+1 : k*in1) = ...
%                         NaN(max(max(K1),max(K2)),in1);
%             end
%         end % Non-reducing end loop
        % Loop for Reducing end angle fitting
        for l = 2%K2
            if any(Mu2{l})
                if size(B,2) < 4
                    for ii = 1 : l
                        tmp = triu(ones(size(B,2)));
                        tmp(tmp == 1) = ...
                            [0 L2{l}(ii,1) 0 L2{l}(ii,2) L2{l}(ii,3) 0];
                        LAM(:,:,ii) = tmp + tril(tmp',-1);
                    end
                elseif size(B,2) < 3
                    LAM = L2{l};                   
                else
                    error('DimensionNotMatch');
                end
                 S_r = struct('Mu',...
                          vmm_ang2rad(reshape(reshape(Mu2{l}(r,:),...
                          [1,size(B,2),l]),[size(B,2) l])'),'Kappa',...
                          Ka2{l},'Lambda',LAM);
                tic
                vmm_r = fitvmmdist(B,l,'Options',options,'Start',...
                            S_r,'Cortype','Sine');
                toc
            else
                error('No convergent initial guesses in disaccharides');
            end
            % Data collection
            Param{r}(1, in1*max(K1)+(l-1)*in2+1 : in1*max(K1)+l*in2) = ...
                vmm_r.Ncomponents;
            if vmm_r.Converged
                Table{r}(:,max(K1)+l) = [vmm_r.Ncomponents vmm_r.AIC ...
                                   vmm_r.BIC vmm_r.Converged vmm_r.Iters]';
                if size(B,2) < 4
                    Param{r}(2 : size(vmm_r.Mu,1) + 1, ...
                        in1*max(K1) + (l-1)*in2 + 1 : ...
                        in1*max(K1) + (l-1)*in2 + 7) = ...
                        [vmm_rad2ang(vmm_r.Mu) vmm_r.Pcomponents' vmm_r.Kappa];
                    for jj = 1 : l
                        Param{r}(size(vmm_r.Mu,1) + jj+1, ...
                        in1*max(K1)+(l-1)*in2 + 8 : in1*max(K1) + l*in2) = ...
                        nonzeros(triu(vmm_r.Lambda(:,:,l)))';
                    end
                elseif size(B,2) < 3
                    Param{r}(2 : size(vmm_r.Mu,1)+1, ...
                        in1*max(K1) + (l-1)*in2+1 : in1*max(K1) + l*in2) = ...
                        [vmm_rad2ang(vmm_r.Mu) vmm_r.Pcomponents' ...
                        vmm_r.Kappa vmm_r.Lambda];
                end
             else
                Table{r}(:,max(K1)+l) = NaN(5,1);
                Param{r}(2:max(max(K1),max(K2))+1, in1*max(K1) + ...
                        (l-1)*in2+1 : in1*max(K1) + l*in2) = ...
                        NaN(max(max(K1),max(K2)),in2);
             end % Convergency 
        end % Reducing End Loop
    end % Repeats
    %matlabpool close;
    
    % Write to files
    pfile = strcat(dir,name,'_','P_',num2str(N),'_',num2str(sd(j)),'.txt');
    tfile = strcat(dir,name,'_','T_',num2str(N),'_',num2str(sd(j)),'.txt');            
    for ii = 1 : N
        dlmwrite(tfile,Table{ii},'-append','delimiter','\t');
        dlmwrite(pfile,Param{ii},'-append','delimiter','\t');
    end
end % Seeds
end % Function

