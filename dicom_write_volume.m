function dicom_write_volume(Volume,filename,volscale,info, startSliceNum)
% This function DICOM_WRITE_VOLUME will write a Matlab 3D volume as
% a stack of 2D slices in separate dicom files.
%
% dicom_write_volume(Volume,Filename,Scales,Info)
%
% inputs,
%   Volume: The 3D Matlab volume
%   Filename: The name of the dicom files
%   Scales: The dimensions of every voxel/pixel
%   Info: A struct with dicom tags and values
%
% Function is written by D.Kroon University of Twente (May 2009)
% modified by Anton Nikolaev RadboudUMC ( Aug 2018)

% Check inputs
if(exist('filename','var')==0), filename=[]; end
if(exist('info','var')==0), info=[]; end
if(exist('volscale','var')==0), volscale=[1 1 1]; end

% Show file dialog if no file name specified
if(isempty(filename))
    [filename, pathname] = uiputfile('*.dcm', 'Save to Dicom');
    filename= [pathname filename];
end

% Add dicom tags to info structure
if(~isstruct(info))
    info=struct;
    % Make random series number
    SN=round(rand(1)*1000);
    % Get date of today
    today=[datestr(now,'yyyy') datestr(now,'mm') datestr(now,'dd')];
    info.SeriesNumber=SN;
    info.AcquisitionNumber=SN;
    info.StudyDate=today;
    info.StudyID=num2str(SN);
    info.PatientID=num2str(SN);
    info.PatientPosition='HFS';
    info.AccessionNumber=num2str(SN);
    info.StudyDescription=['StudyMAT' num2str(SN)];
    info.SeriesDescription=['StudyMAT' num2str(SN)];
    info.Manufacturer='Matlab Convert';
    info.SliceThickness=volscale(3);
    info.PixelSpacing=volscale(1:2);
    info.SliceLocation=(startSliceNum-1)*volscale(3);% modified
end

% Remove filename extention
pl=find(filename=='.'); if(~isempty(pl)), filename=filename(1:pl-1); end

% Read Volume data
disp('Writing Dicom Files...');
for slicenum=startSliceNum:(size(Volume,3)+startSliceNum)
    filenamedicom=[filename number2string(slicenum) '.dcm'];
    % Add slice specific dicom info
    info.InstanceNumber = slicenum;
    info.SliceLocation = info.SliceLocation+volscale(3);
    % Write the dicom file
    %disp(['Writing : ' filenamedicom]);
    %dicomwrite(Volume(:,:,slicenum), colormap(gray), filenamedicom) ;
    dicomwrite(Volume(:,:,slicenum), filenamedicom, info) ;
	if(slicenum==1)
		info2=dicominfo(filenamedicom);
        info.StudyInstanceUID=info2.StudyInstanceUID;
        info.SeriesInstanceUID=info2.SeriesInstanceUID;
	end
end


function numstr=number2string(num)
    num=num2str(num);
    numzeros='000000';  
    numstr=[numzeros(length(num):end) num];