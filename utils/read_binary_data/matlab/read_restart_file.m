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
filenamei = input("Name of the binary restart file written by CaNS [fld.bin]: ");
if isempty(filenamei)
    filenamei = "fld.bin";
end
data       = zeros([ng(1),ng(2),ng(3),4]); % u(:,:,:),v(:,:,:),w(:,:,:),p(:,:,:)
fldinfo    = zeros([2,1]);
[filepath,name,ext] = fileparts(filenamei);
fieldnames = {'u','v','w','p'};
prefix = fullfile(filepath,name);
for q = 1:4
    suffix = ['_',fieldnames{q}];
    if endsWith(name,suffix)
        prefix = fullfile(filepath,name(1:end-length(suffix)));
    end
end
split_files = cell(1,4);
is_split = true;
for q = 1:4
    split_files{q} = [prefix,'_',fieldnames{q},'.bin'];
    is_split = is_split && exist(split_files{q},'file');
end
if is_split
    for q = 1:4
        f = fopen(split_files{q});
        data(:,:,:,q) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
        if q == 1
            fldinfo(:) = fread(f,2,precision);
        end
        fclose(f);
    end
else
    f = fopen(filenamei);
    for q = 1:4
        data(:,:,:,q) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
    end
    fldinfo(:) = fread(f,2,precision);
    fclose(f);
end
%
% store data in arrays
%
u = data(:,:,:,1);
v = data(:,:,:,2);
w = data(:,:,:,3);
p = data(:,:,:,4);
time  =       fldinfo(1);
istep = round(fldinfo(2));
