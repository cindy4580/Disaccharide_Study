clc
clear
mu = zeros(1,2)';
mu(1) = 0;
sigma = ones(1,2)';
p = [ 0.75 0.25];

for j = 1 :4
mu(2) = j;
obj = gmdistribution(mu,sigma,p);
end

