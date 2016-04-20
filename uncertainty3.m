function [ N, total ] = uncertainty3(file,r,S,num,xpos,ypos,flag,pos,pat)
% filtered subsets from trisaccharide MC simulations under density threshold;
% Match local minima into the geometrical positions of circle or ellipse
% Consider local minima from different size of samples
% INPUT VARIABLES:
%       file     - Subsets from density matrix without according directory
%       ub,lb    - Threshold used to get the subset: ub - upper bound; lb - lower bound
%       r        - Energy box size
%       S        - Seed number on the base of 10^4;
%       num      - How many samples to test
%		x/ypos	 - Text box position
%       flag     - Data without omega: 1 ; otherwise: 0
%       pat      - reducing end pattern : 0 or 1 
% Copy-left: Cindy Lee,2014

%% Initialization
str    = regexp(file,'\.','split'); name   = char(str(1));              % Trisaccharide Name
dir    = strcat('~/Data/trimer/',name,'/');                             % File Directory
txt    = ['Sample size: ',num2str(S/10),' * 10^5'];
TT     = 'Local minima prediction';
mini  = [dir,name,'.',num2str(r),'.',num2str(S),'.uncertainty.',pos,'.eps'];
minib = [dir,name,'.',num2str(r),'.',num2str(S),'.uncertainty.',pos,'b.eps'];
M      = load(strcat(dir,file));
temp   = [ ];
%% Data Process
switch pos
    case 'n'
        x = M(1:10^5,2);	y = M(1:10^5,3);	
        if flag
            % Loading Data
            switch S
                case 10
                    for i = 1 : num
                        load(strcat(dir,'/energy/10s.',num2str(r),'.',num2str(i),'_10','.n.mini.mat'));
                        data10((data10(:,6) > 4.11),:) = [ ];
                        temp = [temp; data10];
                    end % for i
                case 20
                    for i = 1 : num
                        load(strcat(dir,'/energy/20s.',num2str(r),'.',num2str(i),'_20','.n.mini.mat'));
                        data20((data20(:,6) > 4.53),:) = [ ];
                        temp = [temp; data20];
                    end % for i
                case 30
                    for i = 1 : num
                        load(strcat(dir,'/energy/30s.',num2str(r),'.',num2str(i),'_30','.n.mini.mat'));
                        data30((data30(:,6) > 4.77),:) = [ ];
                        temp = [temp; data30];
                    end % for i
                case 50
                    for i = 1 : num
                        load(strcat(dir,'/energy/50s.',num2str(r),'.',num2str(i),'_50','.n.mini.mat'));
                        data50((data50(:,6) > 5.07),:) = [ ];
                        temp = [temp; data50];
                    end % for i
                case 100
                    for i = 1 : num
                        load(strcat(dir,'/energy/100s.',num2str(r),'.',num2str(i),'_100','.n.mini.mat'));
                        data100((data100(:,6) > 5.49),:) = [ ];
                        temp = [temp; data100];
                    end % for i
            end % switch S
            % Local minia counts
            [RV, NR, ~, ~] = repval(temp(:,4:5));
            U  = unique(temp(:,4:5),'rows'); % Unique matrix
            U1 = setdiff(U,RV,'rows');
            N  = [ U1 ones(size(U1,1),1); RV NR ];
            total = size(N,1);
            % Visualization
            % background plot
            switch S
                case 10
                    scatter(x,y,0.5,'markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8]);	box on;	hold on;
                case 20
                    scatter(x,y,0.5,'markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7]);   box on; hold on;
                case 30	
                    scatter(x,y,0.5,'markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6]);   box on; hold on;
                case 50	
                    scatter(x,y,0.5,'markerfacecolor',[0.45 0.45 0.45],'markeredgecolor',[0.45 0.45 0.45]); box on; hold on;
            end %switch S
            % hits plot
            for j = 1 : total
                if N(j,3) == 1 || N(j,3) == 2 
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','y');
                    daspect([1,1,1]);   hold on;
                elseif N(j,3) >= 5 && N(j,3) <= 7
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','b');
                    daspect([1,1,1]);   hold on;
                elseif N(j,3) >= 8
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','k');
                    daspect([1,1,1]);   hold on;
                else
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','r');
                    daspect([1,1,1]);   hold on;
                end %if
            end % for
            xlim([0 360]);	ylim([0 360]);
            text(xpos,ypos,txt,'horiz','center','FontSize',14);
            xlabel('\phi(degree)','FontSize',14);	ylabel('\psi(degree)','FontSize',14);
            title(TT,'FontSize',14);
            print(mini,'-depsc2');
        else
            z = M(1:10^5,4);
            zpos = 180;
            % Loading Data
            switch S
                case 10
                    for i = 1 : num
                        load(strcat(dir,'/energy/10s.',num2str(r),'.',num2str(i),'_10','.n.mini.mat'));
                        data10((data10(:,7) > 4.11),:) = [ ]; % more than 20 hits in a box
                        temp = [temp; data10];
                    end % for i
                case 20
                    for i = 1 : num
                        load(strcat(dir,'/energy/20s.',num2str(r),'.',num2str(i),'_20','.n.mini.mat'));
                        data20((data20(:,7) > 4.53),:) = [ ]; % more than 20 hits in a box
                        temp = [temp; data20];
                    end % for i
                case 30
                    for i = 1 : num
                        load(strcat(dir,'/energy/30s.',num2str(r),'.',num2str(i),'_30','.n.mini.mat'));
                        data30((data30(:,7) > 4.77 ),:) = [ ]; % more than 20 hits in a box
                        temp = [temp; data30];
                    end % for i
                case 50
                    for i = 1 : num
                        load(strcat(dir,'/energy/50s.',num2str(r),'.',num2str(i),'_50','.n.mini.mat'));
                        data50((data50(:,7) > 5.07),:) = [ ];
                        temp = [temp; data50];
                    end % for i
                case 100
                    for i = 1 : num
                        load(strcat(dir,'/energy/100s.',num2str(r),'.',num2str(i),'_100','.n.mini.mat'));
                        data100((data100(:,7) > 6.45),:) = [ ];
                        temp = [temp; data100];
                    end % for i
            end % switch
            % Local minia counts
            [RV, NR, ~, ~] = repval(temp(:,4:6));
            U  = unique(temp(:,4:6),'rows'); % Unique matrix
            U1 = setdiff(U,RV,'rows');
            N  = [ U1 ones(size(U1,1),1); RV NR ];
            total = size(N,1);
            % Visualization
            % background plot
            switch S
                case 10
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8]);    box on;	axis equal;
                case 20
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7]);    box on; axis equal;
                case 30	
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6]);    box on;	axis equal;
                case 50	
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.45 0.45 0.45],'markeredgecolor',[0.45 0.45 0.45]);  box on;	axis equal;
            end %switch
            xlim([0 360]);	ylim([0 360]);	zlim([0 360]);	
            text(xpos,ypos,zpos,txt,'horiz','center','FontSize',14);
            xlabel('\phi(degree)','FontSize',14);	ylabel('\psi(degree)','FontSize',14);   
            zlabel('\omega(degree)','FontSize',14);
            title(TT,'FontSize',14);
            print(minib,'-depsc2');
            close all;
            % hits plot
            for j = 1 : total
                if N(j,4) == 1 || N(j,4) == 2 
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[1 1 0]); daspect([1,1,1]);	box on;	
                elseif N(j,4) >= 5 && N(j,4) <= 7
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[0 0 1]); daspect([1,1,1]);	box on;
                elseif N(j,4) >= 8
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[0 0 0]); daspect([1,1,1]);	box on;
                else
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[1 0 0]); daspect([1,1,1]);	box on;
                end %if
            end % for
            xlim([0 360]);	ylim([0 360]);	zlim([0 360]);
            text(xpos,ypos,zpos,txt,'horiz','center','FontSize',14);
            xlabel('\phi(degree)','FontSize',14);	ylabel('\psi(degree)','FontSize',14);   zlabel('\omega(degree)','FontSize',14);
            title(TT,'FontSize',14);
            print(mini,'-depsc2');
        end %if flag
    case 'r'
        switch pat
            case 1
                x = M(1:10^5,4);	y = M(1:10^5,5);
            case 0
                x = M(1:10^5,5);	y = M(1:10^5,6);
        end % switch pat
        if flag
            % Loading Data
            switch S
                case 10
                    for i = 1 : num
                        load(strcat(dir,'/energy/10s.',num2str(r),'.',num2str(i),'_10','.r.mini.mat'));
                        data10((data10(:,6) > 4.11),:) = [ ];
                        temp = [temp; data10];
                    end % for i
                case 20
                    for i = 1 : num
                        load(strcat(dir,'/energy/20s.',num2str(r),'.',num2str(i),'_20','.r.mini.mat'));
                        data20((data20(:,6) > 4.53),:) = [ ];
                        temp = [temp; data20];
                    end % for i
                case 30
                    for i = 1 : num
                        load(strcat(dir,'/energy/30s.',num2str(r),'.',num2str(i),'_30','.r.mini.mat'));
                        data30((data30(:,6) > 4.77),:) = [ ];
                        temp = [temp; data30];
                    end % for i
                case 50
                    for i = 1 : num
                        load(strcat(dir,'/energy/50s.',num2str(r),'.',num2str(i),'_50','.r.mini.mat'));
                        data50((data50(:,6) > 5.07),:) = [ ];
                        temp = [temp; data50];
                    end % for i
                case 100
                    for i = 1 : num
                        load(strcat(dir,'/energy/100s.',num2str(r),'.',num2str(i),'_100','.r.mini.mat'));
                        data100((data100(:,6) > 5.49),:) = [ ];
                        temp = [temp; data100];
                    end % for i
            end % switch S
            % Local minia counts
            [RV, NR, ~, ~] = repval(temp(:,4:5));
            U  = unique(temp(:,4:5),'rows'); % Unique matrix
            U1 = setdiff(U,RV,'rows');
            N  = [ U1 ones(size(U1,1),1); RV NR ];
            total = size(N,1);
            % Visualization
            % background plot
            switch S
                case 10
                    scatter(x,y,0.5,'markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8]);	box on;	hold on;
                case 20
                    scatter(x,y,0.5,'markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7]);   box on; hold on;
                case 30	
                    scatter(x,y,0.5,'markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6]);   box on; hold on;
                case 50	
                    scatter(x,y,0.5,'markerfacecolor',[0.45 0.45 0.45],'markeredgecolor',[0.45 0.45 0.45]); box on; hold on;
            end %switch S         
            % hits plot
            for j = 1 : total
                if N(j,3) == 1 || N(j,3) == 2 
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','y');
                    daspect([1,1,1]);   hold on;
                elseif N(j,3) >= 5 && N(j,3) <= 7
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','b');
                    daspect([1,1,1]);   hold on;
                elseif N(j,3) >= 8
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','k');
                    daspect([1,1,1]);   hold on;
                else
                    rectangle('Position',[N(j,1), N(j,2), r, r],'linestyle','none','FaceColor','r');
                    daspect([1,1,1]);   hold on;
                end %if
            end % for
            xlim([0 360]);	ylim([0 360]);
            text(xpos,ypos,txt,'horiz','center','FontSize',14);
            xlabel('\phi(degree)','FontSize',14);	ylabel('\psi(degree)','FontSize',14);
            title(TT,'FontSize',14);
            print(mini,'-depsc2');
        else
            switch pat
            case 1
                z = M(1:10^5,6);
            case 0
                z = M(1:10^5,7);
            end % switch pat
                      zpos = 180;
            % Loading Data
            switch S
                case 10
                    for i = 1 : num
                        load(strcat(dir,'/energy/10s.',num2str(r),'.',num2str(i),'_10','.r.mini.mat'));
                        data10((data10(:,7) > 4.11),:) = [ ]; % more than 20 hits in a box
                        temp = [temp; data10];
                    end % for i
                case 20
                    for i = 1 : num
                        load(strcat(dir,'/energy/20s.',num2str(r),'.',num2str(i),'_20','.r.mini.mat'));
                        data20((data20(:,7) > 4.53),:) = [ ]; % more than 20 hits in a box
                        temp = [temp; data20];
                    end % for i
                case 30
                    for i = 1 : num
                        load(strcat(dir,'/energy/30s.',num2str(r),'.',num2str(i),'_30','.r.mini.mat'));
                        data30((data30(:,7) > 4.77 ),:) = [ ]; % more than 20 hits in a box
                        temp = [temp; data30];
                    end % for i
                case 50
                    for i = 1 : num
                        load(strcat(dir,'/energy/50s.',num2str(r),'.',num2str(i),'_50','.mini.mat'));
                        data50((data50(:,7) > 5.07),:) = [ ];
                        temp = [temp; data50];
                    end % for i
                case 100
                    for i = 1 : num
                        load(strcat(dir,'/energy/100s.',num2str(r),'.',num2str(i),'_100','.mini.mat'));
                        data100((data100(:,7) > 6.45),:) = [ ];
                        temp = [temp; data100];
                    end % for i
            end % switch
            % Local minia counts
            [RV, NR, ~, ~] = repval(temp(:,4:6));
            U  = unique(temp(:,4:6),'rows'); % Unique matrix
            U1 = setdiff(U,RV,'rows');
            N  = [ U1 ones(size(U1,1),1); RV NR ];
            total = size(N,1);
            % Visualization
            % background plot
            switch S
                case 10
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.8 0.8 0.8],'markeredgecolor',[0.8 0.8 0.8]);    box on;	axis equal;
                case 20
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7]);    box on; axis equal;
                case 30	
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6]);    box on;	axis equal;
                case 50	
                    scatter3(x,y,z,0.5,'markerfacecolor',[0.45 0.45 0.45],'markeredgecolor',[0.45 0.45 0.45]);  box on;	axis equal;
            end %switch
            xlim([0 360]);	ylim([0 360]);	zlim([0 360]);	
            text(xpos,ypos,zpos,txt,'horiz','center','FontSize',14);
            xlabel('\phi(degree)','FontSize',14);	ylabel('\psi(degree)','FontSize',14);   
            zlabel('\omega(degree)','FontSize',14);
            title(TT,'FontSize',14);
            print(minib,'-depsc2');
            close all;
            % hits plot
            for j = 1 : total
                if N(j,4) == 1 || N(j,4) == 2 
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[1 1 0]); daspect([1,1,1]);	box on;	
                elseif N(j,4) >= 5 && N(j,4) <= 7
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[0 0 1]); daspect([1,1,1]);	box on;
                elseif N(j,4) >= 8
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[0 0 0]); daspect([1,1,1]);	box on;
                else
                    plotcube([r r r],[N(j,1),N(j,2),N(j,3)],1,[1 0 0]); daspect([1,1,1]);	box on;
                end %if
            end % for
            xlim([0 360]);	ylim([0 360]);	zlim([0 360]);
            text(xpos,ypos,zpos,txt,'horiz','center','FontSize',14);
            xlabel('\phi(degree)','FontSize',14);	ylabel('\psi(degree)','FontSize',14);   zlabel('\omega(degree)','FontSize',14);
            title(TT,'FontSize',14);
            print(mini,'-depsc2');
        end %if flag
end %switch
end % function
