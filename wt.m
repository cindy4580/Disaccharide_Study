clc
clear
E = load('./12AA_MM.std.wt10.dat');             % changable
h = 10;                                         % size of box;changable

x = E(:,1);
y = E(:,2);
z = E(:,4);                                     % changable according to different DOF
XX = [min(x):h:max(x)];
YY = [min(y):h:max(y)];
M  = size(XX,2);
N  = size(YY,2);
E  = ones(M,N) * ceil(max(z));

for i = 1 : size(x,1)
    ic = [(x(i)-min(x))/h+1, (y(i)-min(y))/h+1];
    E(ic(1),ic(2)) = z(i);
end

figure(1)
surf(XX,YY,E');
Csize = 500;
pos=round(Csize*(3.65-min(E(:)))/(max(E(:)) - min(E(:))));
%Color1 = flipdim(hot(pos));
Color1 = jet(pos);
Color2 = gray(Csize - pos);
colormap([Color1;Color2]);
colorbar;
xlabel('\phi(degree)','FontName','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf','FontSize',11);
ylabel('\psi(degree)','FontName','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf','FontSize',11);
title('Relative free energy matrix (kcal/mol)','FontName','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf');
view(2)
print -depsc2 2.eps  

