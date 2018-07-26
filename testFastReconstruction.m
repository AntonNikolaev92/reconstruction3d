%{
clear all; close all;
path = 'D:\data\cirs';
% users input
vox = [1,1,1]; % voxel size
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
nl = length(zlin);
layers = cell(nl);
zth = vox(3);
tic
for il = 1:nl
    il
    mask = p(3,:) > (zlin(il) - zth) & p(3,:) < (zlin(il) + zth);
    indx = (1:length(mask)).*mask; indx = indx(indx>0);
    layers{il}.x = p(1,indx);
    layers{il}.y = p(2,indx);
    layers{il}.z = p(3,indx);
    layers{il}.v = amp(indx);
    mask = USGRID.z' > (zlin(il) - zth) & USGRID.z' < (zlin(il) + zth);
    indx = (1:length(mask)).*mask; indx = indx(indx>0);
    layers{il}.xq = USGRID.x(indx);
    layers{il}.yq = USGRID.y(indx);
    layers{il}.zq = USGRID.z(indx);
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
    iz = ceil((layers{il}.zq - zmin)/vox(3));
    for i = 1:length(ix)
        if ix(i)>0 & ix(i)<length(xlin) & iy(i)>0 & iy(i)<length(ylin) & iz(i)>0 & iz(i)<length(zlin)
            USDATA(ix(i), iy(i), iz(i)) =  layers{il}.vq(i);
        end
    end
end
toc
%}

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
        
    
    
