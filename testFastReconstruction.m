
clear all; close all;
path = 'D:\data\cirs';
addpath(genpath('D:\nav\libs\matlab'));
% users input
vox = [0.5,0.5,0.5]; % voxel size

% load all data here
[ p, amp ] = createPointsCloud(path);

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
[~ , indx] = sort(p(3,:));
p = p(:,indx);
amp = amp(indx);
USDATA = reconstructAlongZ( p, amp, xlin, ylin, zlin(100:102), vox(3) );
 
%{

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
    %}    
    
    
