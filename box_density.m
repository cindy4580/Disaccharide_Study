function [ ] = box_density(name,box,dbin)
% Get the box_density histogram for 50 samples
% Get boxplot for each density layer
% Input Variables:
%   name - disaccharide name
%   box  - box size to calculate the BOX density originally
%	dbin - data bin centers in series(e.g. 0.025:0.05: 0.975)
% Output results:
%   boxplot figure

%% Initialization
clc
dir     = strcat('~/Data/dimer/',name,'/density/boxdensity/');				% file location
seed    = load('~/Data/dimer/seed1.dat');									% Generating seeds
n       = size(seed,1);														% Number of samples
N       = size(dbin,2);
num     = (360/box)^2;														% Number of boxes in each sample
C       = zeros(N,50);														% Histogram Matrix
dhmatrix= strcat(dir,name,'.dhmatrix.',num2str(box),'.',num2str(N),'bins.dat');
%% Data manipulation
for i = 1 : n
    file = strcat(dir,name,'.box',num2str(box),'.',num2str(seed(i)),'.dat');
    M    = load(file);
    M1   = reshape(M',1,num);
    M1(M1(:) == 0) = NaN;                                                  % Remove empty boxes
    counts= hist(M1,dbin);
    C(:,i) = counts;
end %i
dlmwrite(dhmatrix,C,'delimiter','\t','newline','Unix');
end %function
