
clear all; close all;
path = 'D:\data\cirs';
addpath(genpath('D:\nav\libs\matlab'));
% users input
vox = [0.5,0.5,0.5]; % voxel size
% load all data here
H = []; % transformation metrixes
USDATA = []; 
USGRID = [];
RIO = [ ];
fullfname = [path,'\','H','.mat']; load(fullfname);
fullfname = [path,'\','USDATA','.mat']; load(fullfname);
fullfname = [path,'\','USGRID','.mat']; load(fullfname);
fullfname = [path,'\','RIO','.mat']; load(fullfname);
if  isempty(USGRID) || isempty(USDATA) || isempty(H), return; end;
ns = USGRID.size(1); nel = USGRID.size(2); nfr = USGRID.size(3);% number of elements, samples and frames
if isempty(RIO), RIO = ones(ns, nel, nfr); end;
% create cloud of points
p0 = [ USGRID.x, USGRID.y, zeros(ns*nel,1), ones(ns*nel,1) ]; p0 = p0';
mask = RIO(:)';
indx = [1:length(mask)].*mask; indx = indx(indx>0);
p = zeros(4,sum(mask)*nfr); amp = zeros(1,sum(mask)*nfr);
for ifr = 1:nfr
    p(:, ((ifr-1)*sum(mask)+1):(ifr*sum(mask))) = H(:,:,ifr)*p0(:,indx);
    I = USDATA(:,:,ifr);
    amptmp = I(:);
    amp(((ifr-1)*sum(mask)+1):(ifr*sum(mask))) = amptmp(indx);
end
% create a grid with set accuracy
xmin = min(p(1,:)); xmax = max(p(1,:)); xlin = xmin:vox(1):xmax;
ymin = min(p(2,:)); ymax = max(p(2,:)); ylin = ymin:vox(2):ymax;
zmin = min(p(3,:)); zmax = max(p(3,:)); zlin = zmin:vox(3):zmax;
[X, Y, Z] = meshgrid(xlin, ylin, zlin);
USGRID.x = X(:); USGRID.y = Y(:); USGRID.z = Z(:);
USGRID.size = [ length(xlin), length(ylin), length(zlin) ];
% generate mask


% use statistical correction

% fast reconstruction 
nl = floor(length(zlin));
layers = cell(nl,1);
zth = vox(3);
tic
ns = p;
[~ , indx] = sort(p(3,:));
p = p(:,indx);
amp = amp(indx);
ip = 1;
np = length(p(3,:));
for il = 1:nl
    il
    while p(3,ip) > ( zlin(il) - zth)
        ip = ip - 1;
        if ip == 0, ip = 1; break; end;
    end
    if ip < np, ip = ip + 1; end;
    layers{il}.x = [];
    layers{il}.y = [];
    layers{il}.z = [];
    layers{il}.v = [];
    ip0 = ip;
    tic
    while (p(3,ip) < (zlin(il) + zth) ) && (p(3,ip) > ( zlin(il) - zth) )
        if (ip - ip0) == length(layers{il}.v)
                layers{il}.x = [ layers{il}.x zeros(1,10000) ];
                layers{il}.y = [ layers{il}.y zeros(1,10000) ];
                layers{il}.z = [ layers{il}.z zeros(1,10000) ];
                layers{il}.v = [ layers{il}.v zeros(1,10000) ];
        end
        layers{il}.x(ip-ip0+1) = p(1,ip);
        layers{il}.y(ip-ip0+1) = p(2,ip);
        layers{il}.z(ip-ip0+1) = p(3,ip);
        layers{il}.v(ip-ip0+1) = amp(ip);
        ip = ip+1;
        if ip > np, break; end;
    end
    layers{il}.x = layers{il}.x(1:(ip-ip0 - 1));
    layers{il}.y = layers{il}.y(1:(ip-ip0 - 1));
    layers{il}.z = layers{il}.z(1:(ip-ip0 - 1));
    layers{il}.v = layers{il}.v(1:(ip-ip0 - 1));
    length(layers{il}.v)
    toc
    layers{il}.xq = USGRID.x;
    layers{il}.yq = USGRID.y;
    layers{il}.zq = ones(length(USGRID.x),1).*zlin(il);
end;
clear('p'); clear('p0'); clear('USDATA');
% perform interpolation
pool = gcp();


parfor il = 1:floor(nl)
    F = scatteredInterpolant(layers{il}.x', layers{il}.y', layers{il}.z', layers{il}.v','linear','none');
    layers{il}.vq = F(layers{il}.xq, layers{il}.yq, layers{il}.zq);
end

USDATA = zeros(length(xlin), length(ylin), length(zlin));
for il = 1:floor(nl)
    ix = ceil((layers{il}.xq - xmin)/vox(1));
    iy = ceil((layers{il}.yq - ymin)/vox(2));
    ip = ceil((layers{il}.zq - zmin)/vox(3));
    for i = 1:length(ix)
        if ix(i)>0 & ix(i)<length(xlin) & iy(i)>0 & iy(i)<length(ylin) & ip(i)>0 & ip(i)<length(zlin)
            USDATA(ix(i), iy(i), ip(i)) =  layers{il}.vq(i);
        end
    end
end
toc

SN=round(rand(1)*1000);
today=[datestr(now,'yyyy') datestr(now,'mm') datestr(now,'dd')];
info.SeriesNumber=SN;
info.AcquisitionNumber=SN;
info.StudyDate=today;
info.StudyID=num2str(SN);
info.PatientID=num2str(SN);
info.PatientPosition='Prone';
info.AccessionNumber=num2str(SN);
info.StudyDescription=['CIRS 073' num2str(SN)];
info.SeriesDescription=['CIRS 073' num2str(SN)];
info.Manufacturer='Radboudumc';
info.SliceThickness=vox(3);
info.PixelSpacing=[vox(1) vox(2)];
info.SliceLocation=0;

path = [path,'\','dicom'];
if exist(path), rmdir(path,'s'); end;
mkdir(path);
fname = 'volumeSlice';
fullfname=[path,'\',fname];
volscale = [1 1 1];
dicom_write_volume(USDATA, fullfname, volscale, info);
        
    
    
