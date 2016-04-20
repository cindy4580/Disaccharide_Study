
name = '43YGS';
pattern = '11';

file = strcat(name,'.std.322.pl.dat');
mat  = strcat(name,'.1.dmatrix.5.dat');
[a,b] = glydat3(file,5,1);
glylinkage3(mat,1,pattern)