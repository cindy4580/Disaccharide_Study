clear;
clc

name = '22MMM';
S    = 1;
DenM = strcat(name,'.s',num2str(S),'.dmatrix.5.dat');
%file= strcat('~/Data/dimer/',name,'/',DenM);
file = strcat('~/Data/trimer/',name,'/',DenM);
M    = load(file);
n    = size(M,1);
w    = size(M,2);
y1   = floor(min(M(:,w)));
y2   = ceil(max(M(:,w)));
p    = [ 20 40 50 ];
c    = [ 'r', 'c', 'g'];
for i = 1 : 3
    ener = zeros(1,p(i)); 
    for j = 1 : p(i)
        size = j * (n / p(i));  
        ener(j) = mean(M(1:size,w));
    end
    %figure(i)
    plot(ener,c(i),'linewidth',1.5); hold on;
    ylim([200 205]);
end  
xlabel('Observations');
ylabel('Potential Energy');
title('Energy covergence during one MC simulation');
legend('Every 500 steps','Every 2500 steps','Every 2000 steps');
