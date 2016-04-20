clear;  clc;
seed = load('~/Data/dimer/seed1.dat');
name = '14BB_YN';
for i = 1 : length(seed)
    glylinkageS(name,1,1,seed(i));
end
