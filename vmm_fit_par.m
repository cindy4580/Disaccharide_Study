function [ ] = vmm_fit_par(name,N,K,M,I,eq,varargin) 
%% von Mises Mixture fitting starter script
%  Initial mean directions are randomly selected within 10 degrees of the
%  inputs; Kappa and Lambda initial guesses are fixed for all components at 
%  10 and 1 respectively
%
%  Inputs:
%           Name - Disaccharide name consistent with the whole system
%			N   -- Repeated fit procedures with different initial guesses
%			K   -- Number of components at most
%           M   -- Different samples from MC simulations
%           I   -- Initial mean directions from energy surface local minima
%           eq  -- Start point after equilibrium
%           V1  -- s: Random number seed to generate random initial guesses
%           V2  -- L: Resampling factor

%  NOTE: This script is modified to enable parallel computing on a multi- 
%		 core cpu with MATLAB 2013a therefore contains outdated built-in
%		 functions

%% Initialization
L   = 1;
s   = 27;
% Parse inputs				
if nargin == 0
	return;
end
if nargin < 6
	error('TooFewInputs');
elseif nargin > 7
	s = varargin{1};
    L = varargin{2};
elseif nargin > 6
    s = varargin{1};
end

if size(I,1) ~= K
    error('ClusterCentersMismatchWithClusterNumber');
end
% Inital parameters
rng(s);
p       = size(I,2);                                       % Data Dimension
Kap     = 10;
Mu      = cell(K,N,p);
tmp_M   = zeros(K*N,p);
MP      = cell(N,K,p);

for i   = 1 : K
    for j   = 1 : p
        Mu{i}(:,j) = randi([I(i,j)-5    I(i,j)+5],N,1);
    end
    tmp_M((i - 1) * N + 1 : i * N,:) = Mu{i};     
end

for j   = 1 : N
    MP{j} = tmp_M(j:N:end,:);     
end

if p == 3
    lam  = [ 0 1 0 1 1 0];
    temp = triu(ones(p));
    temp(temp == 1) = lam;
    lambda = temp + tril(temp',-1);
end

% Directories and results
dir     = strcat('~/Data/dimer/',name,'/');
dir1    = strcat(dir,'density/');
Table   = cell(N,5,K);                              % Algorithm Performance
Param   = cell(N,K+1,(p+1)*(p+2)*K/2);               % Parameters Predicted
%% Data process
seed	= load('~/Data/dimer/seed1.dat');
for j 	= M	
	f	= strcat(dir1,name,'.',num2str(seed(j)),'.d.dat');
	A2M = load(f);
    
    % Select every Lth point in radians [-pi, pi]
	A1	= vmm_ang2rad(A2M(eq:L:end,2:p+1));
    clear A2M;
    matlabpool(3)    
    parfor r = 1 : N  % Repeats loop
    %for r = 1 : N
        Mu_tmp = vmm_ang2rad(MP{r});
        
        for k = 1 : K
            
            if p == 1
                error('1D structure is under development');
            elseif p < 3
                S = struct('Mu', Mu_tmp(1:k,:), 'Kappa',...
                    repmat(repmat(Kap,1,p),k,1), 'Lambda',ones(k,1));
%             elseif p < 4
%                 S = struct('Mu', Mu_tmp(1:k,:), 'Kappa',...
%                     repmat(repmat(Kap,1,p),k,1), 'Lambda',...
%                     repmat(lambda,[1 1 k]));
            else
                error('Higher deminsion is not ready yet');
            end
            
            options = statset('MaxIter',300,'TolFun',1e-5);
            tic
            vmm = fitvmmdist( A1, k,'Options', options,'Start', S, ...
                             'Cortype','Sine'); 
            toc
            
            if vmm.Converged
                Table{r}(:,k) = [vmm.Ncomponents vmm.AIC vmm.BIC ...
                                 vmm.Converged vmm.Iters]';
                Param{r}(1,(p+1)*(p+2)/2 *(k-1)+1 : (p+1)*(p+2)/2*k) = ...
                                 vmm.Ncomponents;
                if p < 3 
                    Param{r}(2 : size(vmm.Mu,1) + 1, ...
                    (p+1)*(p+2)/2 * (k-1) + 1 : (p+1)*(p+2)/2 * k) = ...
                    [ vmm_rad2ang(vmm.Mu) vmm.Pcomponents' ...
                      vmm.Kappa vmm.Lambda ];
                elseif p < 4
                    Param{r}(2 : size(vmm.Mu,1) + 1, ...
                    (p+1)*(p+2)/2*(k-1)+1 : (p+1)*(p+2)/2*(k-1)+2*p+1) = ...
                    [ vmm_rad2ang(vmm.Mu) vmm.Pcomponents' vmm.Kappa ];
                    for jj = 1 : k
                        Param{r}(jj+1, (p+1)*((p+2)/2*(k-1) + 2) : ...
                        (p+1)*(p+2)/2 * k) = ...
                        nonzeros(triu(vmm.Lambda(:,:,jj)))';
                    end
                end
            else
                Table{r}(:,k) = zeros(5,1);
                Param{r}(1:K+1,(p+1)*(p+2)/2 *(k-1)+1:...
                               (p+1)*(p+2)/2*k) = zeros(K+1,(p+1)*(p+2)/2);
            end
            
        end % Components
    end % Repeats
    
    pfile = strcat(dir,name,'_','P_',num2str(N),'_',...
            num2str(seed(j)),'.txt');
    tfile = strcat(dir,name,'_','T_',num2str(N),'_',...
            num2str(seed(j)),'.txt');
   
    for ii = 1 :N
            dlmwrite(tfile,Table{ii},'-append','delimiter','\t');        
            dlmwrite(pfile,Param{ii},'-append','delimiter','\t');
    end    
    matlabpool close;    
end % Seed

