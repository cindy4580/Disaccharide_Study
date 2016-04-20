clear;
clc;
box = [ 5 3];
name = '12AA_MM';
B = zeros(20,50);
for i = 1 : size(box,2)
box_density(name,box(i));
end
for i = 1 : size(box,2)
    A = load(strcat(name,'.dhmatrix.',num2str(box(i)),'.dat'));
    for ii = 1 : 50
        B(:,ii) = B(:,ii)/sum(:,ii);
    end % ii
end % i