function [ USDATA ] = reconstructAlongZ( p, amp, qxlin, qylin, qzlin, zth )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% check if the array is sorted along z-axis
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

while p(3,ip) < ( qzlin(1) - zth)
    ip = ip + 1;
    if ip > np, break; end;
end

for il = 1:nl
    il
    while p(3,ip) > ( qzlin(il) - zth)
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
    length(layers{il}.v)
    toc
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
%for il = 1:nl
%    USDATA(:,:,il) = reshape(layers{il}.vq, nxqlin, nyqlin);
%end

for il = 1:floor(nl)
    ix = ceil((layers{il}.xq - qxlin(1))/vox(1));
    iy = ceil((layers{il}.yq - qylin(1))/vox(2));
    iz = ceil((layers{il}.zq - qzlin(1))/vox(3));
    ip = ceil((layers{il}.zq - qzlin(1))/vox(3));
    for i = 1:length(ix)
        if ix(i)>0 & ix(i)<length(qxlin) & iy(i)>0 & iy(i)<length(qylin) & iz(i)>0 & iz(i)<length(qzlin)
            USDATA(ix(i), iy(i), iz(i)) =  layers{il}.vq(i);
            if ix(i)>0 & ix(i)<length(qxlin) & iy(i)>0 & iy(i)<length(qylin) & ip(i)>0 & ip(i)<length(qzlin)
                USDATA(ix(i), iy(i), ip(i)) =  layers{il}.vq(i);
            end
        end
    end
end
