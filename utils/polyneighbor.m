function nn = polyneighbor(shp, linktype, tollerance)
%POLYNEIGHBOR   Return the neighbors of a polygon.
%  Syntax:
%     nn = polyneighbor(shp, linktype, tollerance)
%  Inputs:
%     shp      --  <npts-by-1 strcut>, polygons;
%     linktype --  <char>, {'rock'}|{'queen'}
%                  [ +  *  + ]      [ *  *  * ]
%                  [ *  o  * ]      [ *  o  * ]
%                  [ +  *  + ]      [ *  *  * ]
%                      rock            queen
%  Output:
%     nn     -- <npts-by-1 cell>, spatial contiguity relationship;
%   Copyright 2012 Tang Jianbo, China.
%   This code may be freely used and distributed, so long as it maintains
%   this copyright line.
%   $Revision: 1.0 $     $Date: 2013/01/01 14:12:20 $

% Checking inputs aarguments;
if nargin<1 || isempty(shp)
    error('Not enough input arguments.');
end
if nargin<2 || isempty(linktype)
    linktype = 'rock';
end
if nargin<3 || isempty(tollerance)
    % minima error that two point can be recognized as different
    tollerance = 0.001;  
end

% Getting spatial neighbors of each polygon
if isstruct(shp)
    if ~isfield(shp, 'BoundingBox')
        error('The first input argument is not valid.');
    end
    if ~isequal(shp(1).Geometry, 'Polygon')
        error('The first input argument is not valid.');
    end
    
    npts = length(shp);
    switch lower(linktype)
        case {'rock',  'r', 4}
            min_point_nums = 2;
        case {'queen', 'q', 8}
            min_point_nums = 1;
        otherwise
            error(['The type ''',num2str(linktype), ''' is not vaild.']);
    end
    
    spw = false(npts,npts);
    for i=1:npts
        for j=(i+1):npts
            poly_A = [[shp(i).X]', [shp(i).Y]'];
            ind = ~isnan(poly_A(:,1)) & ~isnan(poly_A(:,2));
            poly_A = poly_A(ind, :);
            clear ind;
            poly_A = round( poly_A./tollerance );
            poly_B = [[shp(j).X]', [shp(j).Y]'];
            ind = ~isnan(poly_B(:,1)) & ~isnan(poly_B(:,2));
            poly_B = poly_B(ind, :);
            clear ind;
            poly_B = round( poly_B./tollerance );
            
            if ~isempty(poly_A) && ~isempty(poly_B)
                interpts = intersect(poly_A, poly_B, 'rows');
                if size(interpts, 1)>= min_point_nums
                    spw(i,j) = true;
                    spw(j,i) = true;
                end
            end
        end
    end % for
    clear x y poly_A poly_B bnd_pts_nums;
    
    nn =cell(npts,1);
    pointids = (1:npts);
    for i=1:npts
        nn{i} = pointids(spw(i,:));
    end % for
else
    error('The first input argument is not a map structure.');
end

end % polyneighbor


