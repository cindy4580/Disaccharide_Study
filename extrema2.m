function [xymax,smax,xymin,smin] = extrema2(xy,varargin)
% EXTREMA2   Gets the extrema points from a surface
% [XMAX,IMAX,XMIN,IMIN] = EXTREMA2(X) returns extrema elements of the matrix X ignoring NaN's, where
% 	XMAX - maxima points in descending order
%	IMAX - linear indexes of the XMAX
%	XMIN - minima points in descending order
%	IMIN - linear indexes of the XMIN.
% The program uses EXTREMA.
% 
% ALGORITHM  Limitations exist as a mathematical program
% 	The extrema points are searched only through the column, the row and
%	the diagonals crossing each matrix element
%	[XMAX,IMAX,XMIN,IMIN] = EXTREMA2(X,1) does the same but without searching through diagonals
%
% NOTE: To change the linear index to (i,j) use IND2SUB. 
% See also EXTREMA, MAX, MIN

% Written by
% Lic. on Physics Carlos Adri�n Vargas Aguilera
% Physical Oceanography MS candidate
% UNIVERSIDAD DE GUADALAJARA 
% Mexico, 2005
%
% nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish. 
% 2006-11-17 : Accept NaN's.
% 2006-11-22 : Fixed bug in INDX (by JaeKyu Suhr)
% 2007-04-09 : Change name to MAXIMA2, and definition added.
% 2014-08    : Modified by cindy in her own directory for circular data analysis

M = size(xy);
if length(M) ~= 2
	error('Entry must be a matrix.')
end
N = M(2);
M = M(1);

% Search peaks through columns:
[smaxcol,smincol] = extremos(xy);

% Search peaks through rows, on columns with extrema points:
im = unique([smaxcol(:,1);smincol(:,1)]); % Rows with column extrema
[smaxfil,sminfil] = extremos(xy(im,:).');

% Convertion from 2 to 1 index:
smaxcol = sub2ind([M,N],smaxcol(:,1),smaxcol(:,2));
smincol = sub2ind([M,N],smincol(:,1),smincol(:,2));
smaxfil = sub2ind([M,N],im(smaxfil(:,2)),smaxfil(:,1));
sminfil = sub2ind([M,N],im(sminfil(:,2)),sminfil(:,1));

% Peaks in rows and in columns:
smax = intersect(smaxcol,smaxfil);
smin = intersect(smincol,sminfil);

% Search peaks through diagonals?
if nargin==1
% Check peaks on down-up diagonal:
 [iext,jext] = ind2sub([M,N],unique([smax;smin]));
 [sextmax,sextmin] = extremos_diag(iext,jext,xy,1);

 % Check peaks on up-down diagonal:
 smax = intersect(smax,[M; (N*M-M); sextmax]);
 smin = intersect(smin,[M; (N*M-M); sextmin]);

 % Peaks on up-down diagonals:
 [iext,jext] = ind2sub([M,N],unique([smax;smin]));
 [sextmax,sextmin] = extremos_diag(iext,jext,xy,-1);

 % Peaks on columns, rows and diagonals:
 smax = intersect(smax,[1; N*M; sextmax]);
 smin = intersect(smin,[1; N*M; sextmin]);
end

% Extrema points:
xymax = xy(smax);
xymin = xy(smin);

% Descending order:
[temp,inmax] = sort(-xymax); clear temp
xymax = xymax(inmax);
smax = smax(inmax);
[xymin,inmin] = sort(xymin);
smin = smin(inmin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [smax,smin] = extremos(matriz)             % Peaks through columns or rows.
smax = [];  smin = [];
for n = 1:size(matriz,2)
	[temp,imaxfil,temp,iminfil] = extrema(matriz(:,n)); clear temp
	if ~isempty(imaxfil)     % Maxima indexes
		imaxcol = repmat(n,length(imaxfil),1);
		smax = [smax; imaxfil imaxcol];
	end
 	if ~isempty(iminfil)     % Minima indexes
  		imincol = repmat(n,length(iminfil),1);
  		smin = [smin; iminfil imincol];
	end
end


function [sextmax,sextmin] = extremos_diag(iext,jext,xy,A)      % Peaks through diagonals (down-up A=-1)
[M,N] = size(xy);
if A==-1
 iext = M-iext+1;
end
[iini,jini] = cruce(iext,jext,1,1);
[iini,jini] = ind2sub([M,N],unique(sub2ind([M,N],iini,jini)));
[ifin,jfin] = cruce(iini,jini,M,N);
sextmax = [];
sextmin = [];
for n = 1:length(iini)
 ises = iini(n):ifin(n);
 jses = jini(n):jfin(n);
 if A==-1
  ises = M-ises+1;
 end
 s = sub2ind([M,N],ises,jses);
 [temp,imax,temp,imin] = extrema(xy(s)); clear temp
 sextmax = [sextmax; s(imax)'];
 sextmin = [sextmin; s(imin)'];
end


function [i,j] = cruce(i0,j0,I,J)
% Indexes where the diagonal of the element io,jo crosses the left/superior
% (I=1,J=1) or right/inferior (I=M,J=N) side of an MxN matrix. 

arriba = 2*(I*J==1)-1;

si = (arriba*(j0-J) > arriba*(i0-I));
i = (I - (J+i0-j0)).*si + J+i0-j0;
j = (I+j0-i0-(J)).*si + J;
