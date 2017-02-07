function [ ] = vmm_fit_plot(name,N,D,K,S,I,L)
%% Visualization of Fitting results
% Based on the structure of results and include the following figures
% AIC/BIC trend
% Scatter plots of fitting results with mean directions
% PDF differences
% Flexible points identifications
%
% Inputs:
%			name -- Disaccharide name code
%			N    -- Repeated fit models with different initial guesses
%			D    -- Resolution of simulation energy/pdf matrix
%           K    -- Largest component number allowed
%			S    -- Monte Carlo seed
%			I	 -- Initial center guesses
%			L    -- Resampling factors

%% Initialization
r    = 32;
c    = [ 0 67 88;
      31 138 112;
      190 219 57;
      225 225 26;
      253 116 0;
      8 217 214;
      225 46 99;
      125 125 125;
      37 42 52]/255;
xpos = 50;
ypos = 100;
seed = load('~/Data/dimer/seed1.dat');
dir  = strcat('~/Data/dimer/',name,'/');
fdir = strcat(dir,'density/');

% Load files
if L == 1
    %f   = strcat(fdir,name,'.std.300.',num2str(seed(S)),'.dat');
    f   = strcat(fdir,name,'.',num2str(seed(S)),'.d.dat');
    p   = strcat(dir,name,'_P_',num2str(N),'_',num2str(seed(S)),'.txt');
    t   = strcat(dir,name,'_T_',num2str(N),'_',num2str(seed(S)),'.txt');
    %d    = strcat(dir,name,'_D_',num2str(N),'_',num2str(seed(S)),'.txt');
else
    f   = strcat(fdir,name,'.',num2str(seed(S)),'.d.dat');
    p   = strcat(dir,name,'_P_',num2str(N),'_',num2str(seed(S)),...
                    '_',num2str(L),'.txt');
    t   = strcat(dir,name,'_T_',num2str(N),'_',num2str(seed(S)),...
                    '_',num2str(L),'.txt');
    %d   = strcat(dir,name,'_D_',num2str(N),'_',num2str(seed(S)),...
    %                '_',num2str(L),'.txt');
end

A2M  = load(f);
P    = load(p);
T    = load(t);
%D1   = load(d);
R    = size(T,1)/N;
scrsz   = get(groot,'ScreenSize');
%% AIC/BIC trends
hp1=figure(1);
for i = 1 : N
    if ~all(T((i-1)*R + 4,:))
        T((i-1)*R+2 : i*R-2, T((i-1)*R+4,:) == 0) = NaN;
    end
    CID = T((i-1)*R + 1,:);
    AIC = T((i-1)*R + 2,:);
    BIC = T((i-1)*R + 3,:);

    plot(CID,AIC,'o-','linewidth',1.5,'markersize',5,...
            'color',[0.5 0.5 0.5]); hold on;
    plot(CID,BIC,'*-','linewidth',1.5,'markersize',5,...
            'color',[0.0 0.0 0.0]);  
end
legend('AIC','BIC');
title('von Mise Mixture model fitting performance','FontSize',14);
xlim([1 K]);  xlabel('Number of fitted clusters');
set(hp1,'Position',[scrsz(3) 1 scrsz(3) scrsz(3)])
ax.XTick = 1 : K;
print(strcat(dir,'FittingPerformance_seed_',num2str(seed(S))),'-dpng');
%% Scatter plot reproduction by von Mises Mixture Model

%%% Lacking consistent colorbar across all subplots %%%

% for ii = 1 : N
% 
%     figure(ii+1)
%     tmp = P((ii - 1) * (K+1) + 1 : ii*(K+1),:);
% 
%     for j = 1 : K 
%         subplot(2,ceil((K)/2),j)
%         if tmp(1,(j-1)*6 + 1) > 0
%             obj = ...
%             vmmdistribution(vmm_ang2rad(tmp(2:j+1,(j-1)*6+1:(j-1)*6+2)),...
%                           tmp(2:j+1,(j-1)*6+4:(j-1)*6+5),tmp(2:j+1,j*6),...
%                           tmp(2:j+1,(j-1)*6 + 3));    
%         
%             z = pdf(vmm_ang2rad(A2M(:,2:3)),obj);
%             z = z/max(z);
%             colormap jet;
%             fastscatter(A2M(:,2), A2M(:,3), z); 
%             xlabel('\phi'); ylabel('\psi'); 
%             xlim([0 360]);  ylim([0 360]); box on;
%         else
%             text(xpos-40,ypos,'not converge');
%             text(xpos-40,ypos+50,strcat(['For ' num2str(j)],...
%                '-cluster model,'));
%             xlim([0 360]);  ylim([0 360]); box on; 
%         end
%         
%         
%     end
%     suptitle('Simulation Pattern Reproduced by von Mises Mixture Model');
%     %saveas
% end
%% Scatter plot with pdf contour
for iii = 1 : N
    hp2 = figure(iii + N +1);
    tmp = P((iii - 1) * (K+1) + 1 : iii*(K+1),:);
    [XX,YY] = meshgrid(0:D:(360-D),0:D:(360-D));
    Temp = vmm_ang2rad([XX(:) YY(:)]);
    for j = 1 : K 
        subplot(2,ceil(K/2),j)
        
        if tmp(1,(j-1)*6 + 1) > 0
            obj = ...
            vmmdistribution(vmm_ang2rad(tmp(2:j+1,(j-1)*6+1:(j-1)*6+2)),...
                          tmp(2:j+1,(j-1)*6+4:(j-1)*6+5),tmp(2:j+1,j*6),...
                          tmp(2:j+1,(j-1)*6 + 3));
        
            z = pdf(Temp,obj)/(D^2);
            z = z/sum(z(:));
            z = reshape(z,360/D,360/D);
            scatter(A2M(1:50:end,2),A2M(1:50:end,3),0.5,...
                      'markerfacecolor',[0.6 0.6 0.6],...
                      'markeredgecolor',[0.6 0.6 0.6]);   
            hold on;    colormap jet;
            contour(XX,YY,z,15); 
            scatter(tmp(2:j+1,(j-1)*6+1),tmp(2:j+1,(j-1)*6+2),...
                          'k+','linewidth',1.5);
            text(xpos,ypos,strcat(num2str(j),' components'),...
                          'Fontsize',12);
            xlim([0 360]);  ylim([0 360]); box on;
            xlabel('\phi','Fontsize',14); ylabel('\psi','FontSize',14);
        else
            text(xpos-40,ypos,'not converge');
            text(xpos-40,ypos+50,strcat(['For ' num2str(j)],...
               '-cluster model,'));
            xlim([0 360]);  ylim([0 360]); box on; 
        end
    end
    
    h = suptitle(strcat(name,': von Mises Mixture Model contour'));
    set(h,'interpreter','none');
    set(hp2,'Position',[scrsz(3) 1 scrsz(3) scrsz(3)])
    print(strcat(dir,'vMMcontour_',num2str(seed(S)),'_',num2str(iii)),...
        '-dpng');
end % Repeats
%% Scatter plot for component centers
m = ['*','o','^','s','p','d','+'];

hp3=figure(2*N+2);
for l = 1 : N
    tmp = P((l - 1) * (K+1) + 1 : l*(K+1),:);
        for ll = 1 : K         
            if tmp(1,(ll-1)*6 + 1) > 0
                scatter(tmp(2:ll+1,(ll-1)*6+1),tmp(2:ll+1,(ll-1)*6+2),...
                'marker',m(ll),'MarkerEdgeColor',c(l,:),...
                'MarkerfaceColor',c(l,:)); hold on;
            end
        end
        box on;
end
h = title(strcat(name,': von Mises Mixture Model centers'));
set(h,'interpreter','none');
x = xlim;   y = ylim;
if x(:,2) > 360
    x(:,2) = 360;
    xlim(x);
end

if y(:,2) > 360
    y(:,2) = 360;
    ylim(y);
end
set(hp3,'Position',[scrsz(3) 1 scrsz(3) scrsz(3)])
print(strcat(dir,'CenterScatter_',num2str(seed(S))),'-dpng');

hp4=figure(2*N + 3);

for ll = 1 : K
    subplot(2,ceil(K/2),ll)
    for l = 1 : N
        tmp = P((l - 1) * (K+1) + 1 : l*(K+1),:);
            scatter(tmp(2:ll+1,(ll-1)*6+1),tmp(2:ll+1,(ll-1)*6+2),...
                'marker',m(ll),'MarkerEdgeColor',c(l,:),...
                'MarkerfaceColor',c(l,:)); hold on;
            ylim(y);
            xlim(x);   
            box on;
    end
end
h = suptitle(strcat(name,': von Mises Mixture Model Centers'));
set(h,'interpreter','none');
set(hp4,'Position',[scrsz(3) 1 scrsz(3) scrsz(3)])
print(strcat(dir,'CenterScatterVsCluster_',num2str(seed(S))),'-dpng');
hold off;

%% Stacked bar plot for mixing weights
%nicer plot would be done in Excel
temp = zeros(size(P,1), K);
for i = 1 : K 
    temp(:,i) = P(:, (i - 1) * 6 + 3);
end
temp(1: K + 1 : end) = [ ];
A   = reshape(temp,[],K);
%fn  = strcat(name,'_WStack_',num2str(seed(S)),'.txt');

hp5=figure(2*N+4);
for k = 2 : K
    B = reshape(A(:,k),K,[]);
    B = B(1:k,:)';
    %dlmwrite(fn,B,'-append','delimiter','\t');
    subplot(2,ceil((K-1)/2),k-1)
    bar(B,'stacked')
    title(strcat(num2str(k),' components'));
    ylim([0 1.1])
end
h = suptitle(strcat(name,': von Mises Mixture Model Weights'));
set(h,'interpreter','none');
set(hp5,'Position',[scrsz(3) 1 scrsz(3) scrsz(3)])
print(strcat(dir,'StackedWeights_',num2str(seed(S))),'-dpng');
%% Initial guesses and fitting results distance
% Here we use euclidean distance first
rng(r);
Mu = zeros(N,2*K);
for i = 1 : K
    Mu(:,(i-1)*2 +1 : (i - 1) * 2 + 2) = ...
        [randi([I(i ,1) - 5 I(i,1) + 5],N,1) ...
        randi([I(i ,2) - 5 I(i,2) + 5],N,1)];
end    
MD = zeros(sum(1:K),N);                                 % Distance matrix  
for i =  1 : N
    hp6 = figure(2*N+ 5 + i);
    for l = 1 : K
        subplot(2,ceil(K/2),l)
        tmp = [Mu(i,1:2:end)' Mu(i,2:2:end)'];
        F = P((i-1) * (K+1) + 2 : (i-1) * (K+1) + (l + 1),  ...
            (l - 1) * 6 + 1 : (l - 1)*6 + 2);
        if ~all(F)
            text(xpos-40,ypos,'not converge');
            text(xpos-40,ypos+50,strcat(['For ' num2str(l)],...
               '-cluster model,'));
            xlim([0 360]);  ylim([0 360]); box on; 

        else
            MD((l^2 -l+2)/2 : (l^2 -l+2)/2 + l-1,i) = ...
                              diag(pdist2(tmp(1:l,:),F));
            line([tmp(1:l,1) F(1:l,1)]',[tmp(1:l,2) F(1:l,2)]',...
                  'color','b','linewidth',1.5); 
            hold on;
            scatter(F(:,1),F(:,2),[],c(1:l,:),'<','filled')
            xlim([0 360]);   ylim([0 360]); box on;
        end
    end
    hold off;
    set(hp6,'Position',[scrsz(3) 1 scrsz(3) scrsz(3)])
    print(strcat(dir,'Travelguess_',num2str(seed(S)),'_',num2str(i)),...
        '-dpng');
    h = suptitle(strcat(name,': Distance from Initial guesses: rep# ',...
                        num2str(i)));
    set(h,'interpreter','none');
end
%% Center and Weights indicator plot
mm = ['k','b','g','y','m','r'];
for i = 1 : N
    hp7 = figure(3*N + 3 + K + i);
    for l = 1 : K
        subplot(2,ceil(K/2),l)
        tmp = P((i-1) * (K+1) + 2 : (i-1) * (K+1) + (l + 1),  ...
            (l - 1) * 6 + 1 : (l - 1)*6 + 3);
        for j = 1 : l
            if tmp(j,3) < 0.01
                scatter(tmp(j,1),tmp(j,2),mm(1),'filled');
            elseif tmp(j,3) < 0.05
                scatter(tmp(j,1),tmp(j,2),mm(2),'filled');
            elseif tmp(j,3) < 0.1
                scatter(tmp(j,1),tmp(j,2),mm(3),'filled');
            elseif tmp(j,3) < 0.3
                scatter(tmp(j,1),tmp(j,2),mm(4),'filled');
            elseif tmp(j,3) < 0.6
                scatter(tmp(j,1),tmp(j,2),mm(5),'filled');
            elseif tmp(j,3) <= 1
                scatter(tmp(j,1),tmp(j,2),mm(6),'filled');
            end
            hold on;
            ylim(y);xlim(x); box on;
        end
    end
    set(hp7,'Position',[scrsz(3) 1 scrsz(3) scrsz(3)])
    h = suptitle(strcat(name,': von Mises Mixture Model Weights'));
    set(h,'interpreter','none');
    print(strcat(dir,'CenterWeights_',num2str(seed(S)),'_',num2str(i)),...
        '-dpng');
end

%% Differences between pdfs
% for j = 1 : N
%     
%     %figure(j + 2*N + 1)
%     
% end
%% Pdf differences analysis
% Histogram Plot

% for jj = 1 : N
%     figure(jj + 3*N + 1)
%     tmp = D1((jj - 1) * 360/D + 1 : jj*360/D,:);
%     
%     for k = 1 : K - 1
%         D2 = tmp(:,(k-1)*360/D + 1 : k*360/D);
%         subplot(ceil((K-1)/2),2,k)
%         histogram(D2,'Normalization','probability');
%         ID = strcat(num2str(k+1),' components');
%         ylim([0 0.4]);
%         ylabel(ID); box on;
%     end
% 
%     
% end

% Box plot 

% for jjj = 1 :N
%     figure(jjj+ 4*N + 1)
%     tmp = D1((jjj - 1) * 360/D + 1 : jjj*360/D,:);
%     D2  = zeros((360/D)^2,K-1);
%     for k = 1 : K -1
%         D2(:,k) = reshape(tmp(:,(k-1)*360/D + 1 : k*360/D),[(360/D)^2 1]);
%     end
%     boxplot(D2,'Labels',{'K = 2','K = 3','K = 4','K = 5','K = 6'},...
%         'Symbol','.','color',[0.2 0.6 1.0],'OutlierSize',6);
%     h = findobj(gca,'tag','Median');
%     set(h,'color',[1.0 0.4 0.0],'linewidth',1.5);
%     h1 = findobj(gca,'tag','Outliers');
%     set(h1,'markeredgecolor',[1.0 0.2 0.4]);
% %     a = get(get(gca,'children'),'children');
% %     te = get(a,'tag'); 
%     
%     
%     %title('Differences between simulation and predicted pdfs',...
%      %   'FontSize',10)
%       hold on;
%     plot([-1 100],[ 0 0],':k');
%     print('vmmdiffbox','-depsc2');
% end
%% Cluster plot
            % Visualization
            %subplot(2,3,k-1)
%           idx = cluster(A1,vmm);
% 			I(:,(r-1)*(K-1)+(k-1)) = idx;
           	% for m = 1 : k
            %    scatter(A1(idx == m,1),A1(idx == m,2),10,color(m,:)); 
            %   hold on; xlim([-pi pi]); ylim([-pi pi]);
            % end
           	%suptitle('VMM fitting for Man-\alpha(1,2)-man')

end % Function


