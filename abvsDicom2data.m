% abvs cone 2 data
 addpath(genpath('D:\nav\libs\matlab'));
% addpath(genpath('E:\libs\matlab'));
 path = 'D:\data\cirs';
%path = 'E:\data\cirs';
dx = -29.7055; dy = 0; dangle = 45;

R = getRotationMetrix(0,0,dangle);
T = eye(4); T(1,4) = dx; T(2,4) = dy;
Hc = R*T;
% prepare grid
fullfname = [path,'\','DATA#0_0.IMA'];
info = dicominfo(fullfname); info.NumberOfFrames = info.NumberOfFrames - 8;

% USDATA
fullfname = [path, '/', 'DATA_NOISE0','.IMA']; curUSDATA = dicomread(fullfname);
NOISE = uint32(curUSDATA(:,:,1,1));
for ifr = 2:(info.NumberOfFrames - 8)
    NOISE = NOISE + uint32(curUSDATA(:,:,1,ifr));
end
NOISE = uint8( NOISE/(info.NumberOfFrames - 8) );


fnames = dir(path); strfnames = [];
for ifname = 1:length(fnames)
    strfnames = [strfnames,fnames(ifname).name];
end;
fnames = regexp(strfnames,'DATA#\d*_\d*.IMA','match');
strStartAngles = regexp(fnames,'\d*_\d*','match');
for isa = 1:length(strStartAngles)
    startAngles(isa) = str2double( strrep( strStartAngles{isa},'_','.') );
end;

H = zeros(4,4,info.NumberOfFrames*length(startAngles));
USDATA = zeros(info.Height,info.Width, info.NumberOfFrames*length(startAngles));
da = 360/info.NumberOfFrames;
nfr = 0;
for isa = 1:length(startAngles)
    fullfname = [path, '\', fnames{isa}]; curUSDATA = dicomread(fullfname);
    for ia = 1:info.NumberOfFrames
        nfr = nfr + 1;
        R = getRotationMetrix(0,startAngles(isa)+(ia-1)*da ,0);
        H(:,:,nfr) = R*Hc;
        USDATA(:,:,nfr) = curUSDATA(:,:,1,ia) - NOISE + 1;
    end
end

% USGRID
dh = info.PixelSpacing(1); dw = info.PixelSpacing(2);
w = double(info.Width)*dw; h = double(info.Height)*dh;
xlin = linspace(0,w,info.Width); ylin = linspace(0,h,info.Height); zlin = 0;
[ X, Y, Z ] = meshgrid(xlin, ylin, zlin);
USGRID = struct('x',[],'y',[], 'z',[],'size', [1,3]);
USGRID.x = X(:); USGRID.y = Y(:); USGRID.z = Z(:);
USGRID.size = [ uint32(info.Width), uint32(info.Height), uint32(nfr) ];

% RIO
p0 = [USGRID.x, USGRID.y, USGRID.z, ones(length(USGRID.x),1)]; p0 = p0';
p = Hc*p0;
mask = [ p(1,:)>0 ] .* [ p(2,:)>0 ];
RIO = reshape(mask,info.Height,info.Width);

fullfname = [path,'\','USDATA','.mat']; save(fullfname,'-v7.3','USDATA');
fullfname = [path,'\','H','.mat']; save(fullfname,'-v7.3', 'H');
fullfname = [path,'\','USGRID','.mat']; save(fullfname,'-v7.3', 'USGRID');
fullfname = [path,'\','RIO','.mat']; save(fullfname,'-v7.3', 'RIO');



