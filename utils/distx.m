function D = distx(point,pointset)
%DISTX   compute euclidean distance between one point and a point set.
%  Syntax:
%     D = DISTX(point,pointset)
%  Input:
%     point - one observation;
%     pointset - points set;
%  Output:
%     D - distance, the same length as point set;
%
%  See also SQRT
%
%  Copyright 2013, Tang Jianbo, China.
%  This code may be freely used and distributed, so long as it maintains
%  this copyright line.
%  Version: 1.0,   Date: Mar-16-2013 22:30:20
%

if nargin<2
    error('Not enough input arguments.');
end
if size(point,2)~=size(pointset,2)
    error('The two arguments are not the same dimension.');
end
if isempty(point)||isempty(pointset)
    D = [];
    return;
end

[npts,dims] = size(pointset);
if(dims==1)
    D = abs(ones(npts,1)*point-pointset);
elseif(dims>1)
    D = sqrt(sum((ones(npts,1)*point-pointset).^2, 2));
else
    D = [];
end % if
end % function

