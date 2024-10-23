function [s, connect] = redcap_spatialmatrix(shpfile, linktype, tollerance)
if nargin<2 || isempty(linktype)
    linktype = 'rock';
end
if nargin<3 || isempty(tollerance)
    % minima error that two point can be recognized as different
    tollerance = 0.01;  
end

%redcap_spatialmatrix   Return the neighbors of a polygon.
%  Syntax:
%     nn = redcap_spatialmatrix(shpfile, linktype, tollerance)
%  Inputs:
%     shpfile  --  ESRI shapefile;
%     linktype --  <char>, {'rock'}|{'queen'}
%                  [ +  *  + ]      [ *  *  * ]
%                  [ *  o  * ]      [ *  o  * ]
%                  [ +  *  + ]      [ *  *  * ]
%                      rock            queen
%  Output:
%     nn     -- <npts-by-1 cell>, spatial contiguity relationship;

shp = shaperead(shpfile);
fn = fieldnames(shp);
fn = fn(~ismember(fn,{'Geometry','BoundingBox','X','Y'}));

if isempty(fn) || length(shp)<1
    s = [];
    connect = [];
    return;
end

n = length(shp);
m = length(fn);

s.name  = fn;
s.data  = cell(n,m);

for i=1:n
    for j=1:m
        fname = fn{j};
        s.data{i,j} = shp(i).(fname);
    end
end

connect = polyneighbor(shp, linktype, tollerance);

end % redcap_spatialmatrix







