function [ USDATA ] = reconstructAlongZ( p, amp, qxlin, qylin, qzlin, zth )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% check if the array is sorted along z-axis
if isempty(p), return; end;
if isempty(amp), return; end;
if length(qxlin)<2, return; end;
if length(qylin)<2, return; end;
if length(qzlin)<2, return; end;
if zth<0, return; end;

if ~issorted(p(3,:))
    [~ , indx] = sort(p(3,:));
    p = p(:,indx);
    amp = amp(indx);
end
vox(1) = qxlin(2) - qxlin(1);
vox(2) = qylin(2) - qylin(1);
vox(3) = qzlin(2) - qzlin(1);

[X, Y, Z] = meshgrid(qxlin, qylin, qzlin);
USGRID.x = X(:); USGRID.y = Y(:); USGRID.z = Z(:);
USGRID.size = [ length(qxlin), length(qylin), length(qzlin) ];

nl = floor(length(qzlin));
layers = cell(nl,1);
ip = 1;
np = length(p(3,:));

while p(3,ip) < ( qzlin(1) - zth )
    if ip == np, break; end;
    ip = ip + 1;
end

for il = 1:nl
    while p(3,ip) > ( qzlin(il) - zth)
        if ip == 1, break; end;
        ip = ip - 1;
    end
    if ip < np, ip = ip + 1; end;
    layers{il}.x = [];
    layers{il}.y = [];
    layers{il}.z = [];
    layers{il}.v = [];
    ip0 = ip;
    while (p(3,ip) < (qzlin(il) + zth) ) && (p(3,ip) > ( qzlin(il) - zth) )
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
    layers{il}.xq = USGRID.x;
    layers{il}.yq = USGRID.y;
    layers{il}.zq = ones(length(USGRID.x),1).*qzlin(il);
end

% perform interpolation
parfor il = 1:nl
    F = scatteredInterpolant(layers{il}.x', layers{il}.y', layers{il}.z', layers{il}.v','linear','none');
    layers{il}.vq = F(layers{il}.xq, layers{il}.yq, layers{il}.zq);
end

nxqlin = length(qxlin); nyqlin = length(qylin); nzqlin = length(qzlin);
USDATA = zeros(nxqlin, nyqlin, nzqlin);

for il = 1:nl
    ix = floor((layers{il}.xq - qxlin(1))/vox(1))+1;
    iy = floor((layers{il}.yq - qylin(1))/vox(2))+1;
    iz = floor((layers{il}.zq - qzlin(1))/vox(3))+1;
    for i = 1:length(ix)
        if ix(i)>0 & ix(i)<=nxqlin & iy(i)>0 & iy(i)<=nyqlin & iz(i)>0 & iz(i)<=nzqlin
            USDATA(ix(i), iy(i), iz(i)) =  layers{il}.vq(i);
        end
    end
end
