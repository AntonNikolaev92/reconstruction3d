function [ p, amp ] = createPointsCloud( path )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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

end

