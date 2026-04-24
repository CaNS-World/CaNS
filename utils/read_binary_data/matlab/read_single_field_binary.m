% -
%
% SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
% SPDX-License-Identifier: MIT
%
% -
%
% setting up some parameters
%
precision = 'double';     % precision of the real-valued data
r0 = [0.,0.,0.];    % domain origin
%
% read geometry file
%
geofile  = "geometry.out";
data = dlmread(geofile);
ng   = data(1,:);
l    = data(2,:);
dl   = l./ng;
%
% read and generate grid
%
xp = linspace(r0(1)+dl(1)/2.,r0(1)+l(1)-dl(1)/2.,ng(1)); % centered  x grid
yp = linspace(r0(2)+dl(2)/2.,r0(2)+l(2)-dl(2)/2.,ng(2)); % centered  y grid
zp = linspace(r0(3)+dl(3)/2.,r0(3)+l(3)-dl(3)/2.,ng(3)); % centered  z grid
xu = xp + dl(1)/2.; % staggered x grid
yv = yp + dl(2)/2.; % staggered y grid
zw = zp + dl(3)/2.; % staggered z grid
if(exist('grid.bin','file'))
    f   = fopen('grid.bin');
    grid_z = fread(f,[ng(3),4],precision);
    fclose(f);
    zp = r0(3) + grid_z(:,3)'; % centered  z grid
    zw = r0(3) + grid_z(:,4)'; % staggered z grid
end
%
% read checkpoint binary file
%
filenamei = input("Name of the binary file written by CaNS (e.g. vex_fld_0000000.bin)]: ")
if isempty(filenamei)
end
iskipx      = input("Data saved every (ix, iy, iz) points. Value of ix? [1]: ")
if isempty(iskipx)
    iskipx = 1
end
iskipy      = input("Data saved every (ix, iy, iz) points. Value of iy? [1]: ")
if isempty(iskipy)
    iskipy = 1
end
iskipz      = input("Data saved every (ix, iy, iz) points. Value of iz? [1]: ")
if isempty(iskipz)
    iskipz = 1
end
iskip       = [iskipx,iskipy,iskipz]
n           = floor((ng-1)./iskip)+1
f = fopen(filenamei);
fld = fread(f,prod(n),precision);
fclose(f);
if numel(fld) ~= prod(n)
    error("expected %d values for the requested skip, found %d",prod(n),numel(fld))
end
data = reshape(fld,[n(1),n(2),n(3)]);
