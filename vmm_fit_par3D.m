function [ ] = vmm_fit_par3D(name,N,D,K,M,I,varargin) 
%% von Mises Mixture fitting and measuring the pdf difference 
%  Initial mean directions are randomly flexible within 10 degrees; Kappa
%  and Lambda are fixed at 10 and 1 respectively
%
%  Inputs:
%			N   -- Repeated fit models with different initial guesses
%			D   -- Resolution of simulation energy/pdf matrix
%			K   -- Number of components at most
%           M   -- Different samples from MC simulations
%           I   -- Initial mean directions from energy surface local minima
%           V1  -- L: Resampling factor
%           V2  -- s: Random number seed to generate initial guesses
%           V3  -- C: Cutoff for pdf comparison
%
%  NOTE: This script is modified to enable parallel computing on a multi- 
%		 core cpu with MATLAB 2013b 

%% Initialization
R 	= -1.9872041 * 10^(-3);                            % Ideal gas constant
T 	= 300;                                                    % Temperature
C 	= 4.12;											  % >=100 dots in a box
s   = 32; 						   % Random number seed for initial guesses
L   = 1;
% Parse inputs				
if nargin == 0
	return;
end
if nargin < 6
	error('TooFewInputs');
elseif nargin > 8
    C = varargin{3};
elseif nargin > 7
	A2M = varargin{2};
    L = varargin{1};
elseif nargin > 6
    L = varargin{1};
end

if size(I,1) ~= K
    error('ClusterCentersMismatchWithNumber');
end

% [XX,YY] = meshgrid(0:D:(360-D),0:D:(360-D));
% Temp = vmm_ang2rad([XX(:) YY(:)]);

% Inital parameters
rng(s);
p       = size(I,2);
Kap     = 10;
Mu      = cell(K,N,p);
tmp_M   = zeros(K*N,p);
MP      = cell(N,K,p);

for i = 1 : K
    for j = 1 : p
        Mu{i}(:,j) = randi([I(i,j)-5    I(i,j)+5],N,1);
    end
    
    tmp_M((i - 1) * N + 1 : i * N,:) = Mu{i};     
end

for j = 1 : N
    MP{j} = tmp_M(j:N:end,:);     
end

if p == 3
    lam  = [ 0 1 0 1 1 0];
    temp = triu(ones(p));
    temp(temp == 1) = lam;
    lambda = temp + tril(temp',-1);
end

% Results
dir   = strcat('~/Data/dimer/',name,'/');
Table = cell(N,5,K);                                % Algorithm Performance
Param = cell(N,K+1,(p+1)*(p+2)*K/2);                   % Parameters Predicted
%Diff  = cell(N,360/D,360/D * K);                       % Energy difference
%% Data load
%dir1    = strcat('~/Data/dimer/',name,'/density/');
%dir2    = strcat('~/Data/dimer/',name,'/ener/');
seed	= load('~/Data/dimer/seed1.dat');

for j 	= M
	
	%f	= strcat(dir1,name,'.',num2str(seed(j)),'.d.dat');
    %f   = strcat(dir1,name,'.std.300.',num2str(seed(j)),'.dat');
    f    = strcat(dir,name,'.1.dmatrix.5.dat');
	%e   = strcat(dir2,name,'.wt',num2str(D),'.',num2str(seed(j)),'.dat');
	A2M = load(f);
	%Ene = flipud(load(e));
    
    % Select every Lth point in radians [-pi, pi]
	A1	= vmm_ang2rad(A2M(1:1:end,2:p+1));   
	%Ene(Ene(:) == 100) = 0;
	%Ene(Ene(:) == 100) = NaN;
	%prob = exp(Ene/(R*T))/D^2;
    %prob(prob(:) == 1/D^2) = 0;
    
    %matlabpool(3)
    
    %parfor r = 1 : N                                         % Repeats loop
    for r = 1 : N

        Mu_tmp = vmm_ang2rad(MP{r});
        
        for k = 2 % 1: K
            if p < 2
                error('1D structure is under development');
            elseif p < 3
                S = struct('Mu', Mu_tmp(1:k,:), 'Kappa',...
                    repmat(repmat(Kap,1,p),k,1), 'Lambda',ones(k,1));
            elseif p < 4
                S = struct('Mu', Mu_tmp(1:k,:), 'Kappa',...
                    repmat(repmat(Kap,1,p),k,1), 'Lambda',...
                    repmat(lambda,[1 1 k]));
            else
                error('Higher deminsion is not ready yet');
            end
            options = statset('MaxIter',500,'TolFun',1e-5);
            tic
            vmm = fitvmmdist( A1, k,'Options', options,'Start', S, ...
                             'Cortype','Sine'); toc
            
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
                            
                % Predicted pdf if algorithm converges
%                 z = pdf(Temp,vmm);
%                 z = reshape(z,360/D,360/D);12BA_NM_T_9_cmb.txt
%                 d = prob - z/(sum(z(:))*D^2);
% 
%                 Diff{r}(:,(k-1)*360/D + 1 : (k)*360/D) = d;
            else
                Table{r}(:,k) = zeros(5,1);
                Param{r}(1:K+1,(p+1)*(p+2)/2 *(k-1)+1:...
                               (p+1)*(p+2)/2*k) = zeros(K+1,(p+1)*(p+2)/2);
            end
            
        end % Components
    end % Repeats
    
    % Write to files
    if L == 1
        pfile = strcat(dir,name,'_','P_',num2str(N),'_',...
                num2str(seed(j)),'.txt');
        tfile = strcat(dir,name,'_','T_',num2str(N),'_',...
                num2str(seed(j)),'.txt');
        dfile = strcat(dir,name,'_','D_',num2str(N),'_',...
                num2str(seed(j)),'.txt');
    else
       pfile = strcat(dir,name,'_','P_',num2str(N),'_',...
                num2str(seed(j)),num2str(L),'.txt');
       tfile = strcat(dir,name,'_','T_',num2str(N),'_',...
                num2str(seed(j)),num2str(L),'.txt');
       dfile = strcat(dir,name,'_','D_',num2str(N),'_',...
                num2str(seed(j)),num2str(L),'.txt');
    end
    
    for ii = 1 :N
            dlmwrite(tfile,Table{ii},'-append','delimiter','\t');        
            dlmwrite(pfile,Param{ii},'-append','delimiter','\t');
%           dlmwrite(dfile,Diff{ii},'-append','delimiter','\t');
    end
    %matlabpool close;
end % Seed

