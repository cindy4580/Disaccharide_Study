clc
clear
%% Initialization
name1   = '12AA_MM';
name2   = '22MMM';
incr    = 5;
num     = 360/incr;
lname   = 0 : incr : (360 - incr);
ener    = 5;                                        % energy threshold for at least 100 dots in each box
nsam    = 40;                                       % threshold for sample validity
seed1   = load('~/Data/dimer/seed1.dat');
seed2   = load('~/Data/trimer/seed2.dat');
dir1    = strcat('~/Data/dimer/',name1,'/ener/');
dir2    = strcat('~/Data/trimer/',name2,'/ener/');
hm1     = strcat('~/Data/trimer/',name2,'/ener/ditr.dat');
hm1h    = strcat('~/Data/trimer/',name2,'/ener/ditrh.dat');
hm2     = strcat('~/Data/trimer/',name2,'/ener/ditn.dat');
hm2h     = strcat('~/Data/trimer/',name2,'/ener/ditnh.dat');
diener  = cell(1, size(seed1,1));                   % disaccharide energy matrix as reference
trire   = cell(1, size(seed1,1));                   % trisaccharide & reducing end
trine   = cell(1, size(seed1,1));                   % trisaccharide & non-reducing end
disam   = cell(num,num);                            % disaccharide samples
trirsam = cell(num,num);                            % trisaccharide & reducing end sample
trinsam = cell(num,num);                            % trisaccharide & non-reducing end sample
testdr  = cell(num,num);                            % test statistics cell [h,p] : di VS. tri-re
testdrp = NaN(num);                                 % test : P-value matrix
testdrh = NaN(num);                                 % test : decision matrix
testdn  = cell(num,num);                            % test statistics cell [h,p] : di VS. tri-ne
testdnp = NaN(num);                                 % test : P-value matrix
testdnh = NaN(num);                                 % test : decision matrix
cnt_good_di = 0;                                    % # of good boxed ( valid dots > nsam): disaccharide
cnt_good_tr = 0;                                    % # of good boxed : trisaccharide of reducing end
cnt_good_tn = 0;                                    % # of good boxed : trisaccharide of non-reducing end
cnt_nsam_di = 0;                                    % # of never-sampled boxes : disaccharide
cnt_nsam_tr = 0;                                    % # of never-sampled boxes : trisaccharide of reducing end
cnt_nsam_tn = 0;                                    % # of never-sampled boxes : trisaccharide of non-reducing end
%% Data cleaning
% Load and assign NaNs : disaccharide sample
for i = 1 : size(seed1,1)
    file = strcat(dir1,name1,'.wt5.',num2str(seed1(i)),'.dat');
    diener{i} = importdata(file,'\t',0);
    diener{i}(diener{i}(:,:) == 100 ) = NaN; 
end
% Load and assign NaNs : trisaccharide & reducing end
for i   = 1 : size(seed2,1)
    file = strcat(dir2,name2,'.wt5.',num2str(seed2(i)),'.r.dat');
    trire{i} = importdata(file,'\t',0);
    trire{i}(trire{i}(:,:) == 100 ) = NaN; 
end
% Load and assign NaNs : trisaccharide & non-reducing end
for i   = 1 : size(seed2,1)
    file = strcat(dir2,name2,'.wt5.',num2str(seed2(i)),'.n.dat');
    trine{i} = importdata(file,'\t',0);
    trine{i}(trine{i}(:,:) == 100 ) = NaN; 
end
%% diener{i}(:,:) > ener
% collect samples: disaccharide
for i = 1 : num
    for j = 1 : num
        for k = 1 : size(seed1,1)
        disam{i,j}(k) = diener{k}(i,j);
        end
    end
end
% collect samples: trisaccharide & reducing end
for i = 1 : num
    for j = 1 : num
        for k = 1 : size(seed2,1)
        trirsam{i,j}(k) = trire{k}(i,j);
        end
    end
end
% collect samples: trisaccharide & non-reducing end
for i = 1 : num
    for j = 1 : num
        for k = 1 : size(seed2,1)
        trinsam{i,j}(k) = trine{k}(i,j);
        end
    end
end
%% Difference Analysis
% comparison between disaccharide and reducing end linkage of trisaccharide  
for i   = 1 : num
    for j    = 1 : num
        n_sam_di    = size(seed1,1) - sum(isnan(disam{i,j}));
        n_sam_tr    = size(seed2,1) - sum(isnan(trirsam{i,j}));
        n_good_di   = n_sam_di - sum(disam{i,j} > ener);
        n_good_tr   = n_sam_tr - sum(trirsam{i,j} > ener);
        if n_good_di >= nsam && n_good_tr >= nsam
           [h,p] = kstest2(disam{i,j}, trirsam{i,j},'Alpha',0.01);
           testdr{i,j} = [h,p];
           testdrp(i,j) = p;
           testdrh(i,j) = h;
        end
        % Statistics
        if n_good_di >= nsam
            cnt_good_di = cnt_good_di + 1;
        end
        if n_good_tr >= nsam
            cnt_good_tr = cnt_good_tr + 1;
        end
        if n_sam_di == 0
            cnt_nsam_di = cnt_nsam_di + 1;
        end
        if n_sam_tr == 0
            cnt_nsam_tr = cnt_nsam_tr + 1;
        end
    end
end
%HMobj = HeatMap(flipud(testdrp),'RowLabels',lname,'ColumnLabels',lname);
%plot(HMobj);
%print(hm1,'-depsc2');
dlmwrite(hm1,testdrp,'delimiter','\t','newline','Unix');
dlmwrite(hm1h,testdrh,'delimiter','\t','newline','Unix');
pgood_di    = cnt_good_di/(num^2 - cnt_nsam_di);     % pct. of good samples relative to sampled boxes : di
psam_di     = (num^2 - cnt_nsam_di)/num^2;           % pct. of ever-sampled relative to whole boxes : di
pgood_tr    = cnt_good_tr/(num^2 - cnt_nsam_tr);     % pct. of good samples relative to sampled boxes :tri-re
psam_tr     = (num^2 - cnt_nsam_tr)/num^2;           % pct. of ever-sampled relative to whole boxes : tri-re
% comparison between disaccharide and non-reducing end linkage of trisaccharide  
for i   = 1 : num
    for j    = 1 : num
        n_sam_di    = size(seed1,1) - sum(isnan(disam{i,j}));
        n_sam_tn    = size(seed2,1) - sum(isnan(trinsam{i,j}));
        n_good_di   = n_sam_di - sum(disam{i,j} > ener);
        n_good_tn   = n_sam_tn - sum(trinsam{i,j} > ener);
        if n_good_di >= nsam && n_good_tn >= nsam
           [h,p] = kstest2(disam{i,j}, trinsam{i,j},'Alpha',0.01);
           testdn{i,j} = [h,p];
           testdnp(i,j) = p;
           testdnh(i,j) = h;
        end
        if n_good_tn >= nsam
            cnt_good_tn = cnt_good_tn + 1;
        end
        if n_sam_tn == 0
            cnt_nsam_tn = cnt_nsam_tn + 1;
        end
    end
end
%HeatMap(flipud(testdnp),'RowLabels',lname,'ColumnLabels',lname);
%print(hm2,'-depsc2');
dlmwrite(hm2,testdnp,'delimiter','\t','newline','Unix');
dlmwrite(hm2h,testdnh,'delimiter','\t','newline','Unix');
pgood_tn    = cnt_good_tn/(num^2 - cnt_nsam_tn);     % pct. of good samples relative to sampled boxes :tri-re
psam_tn     = (num^2 - cnt_nsam_tn)/num^2;           % pct. of ever-sampled relative to whole boxes : tri-re
