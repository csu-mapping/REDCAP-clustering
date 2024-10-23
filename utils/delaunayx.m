function Egs = delaunayx(x, y, beta, fig)
%   Copyright 2012 Tang Jianbo, China.
%   This code may be freely used and distributed, so long as it maintains
%   this copyright line.
%   $Revision: 1.0 $     $Date: 2013/01/01 14:12:20 $
%DELAUNAYX  
if nargin<3||isempty(beta)
    beta = 2;
end
if nargin<4||isempty(fig)
    fig = 'on';
end
x = x(:);
y = y(:);
npts = length(x);

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);
conner = [[xmin; xmax; xmin; xmax],[ymin; ymin; ymax; ymax]];
x = [x; conner(:,1)];
y = [y; conner(:,2)];
dt = DelaunayTri(x,y);
E = edges(dt);
L = edgesLength(E,x,y);
threshold = mean(L)+std(L)*beta;
Egs = E(L<threshold,:);

connerids = npts+[1,2,3,4];
ind = ismember(Egs(:,1),connerids)|ismember(Egs(:,2),connerids);
Egs = Egs(~ind,:);
switch lower(fig)
    case {'yes', 'on', 'y', 1}
        if ~isempty(Egs)
            figure;
            showedges([x,y],Egs,'m');
            hold on;
            scattx([x,y],[],{'.'},[],{'k'});
            scattx(conner,[],{'+'},[],{'b'});
        end
    otherwise
end
end % function



function varargout = showedges(vertex,edges,varargin)
%SHOWEDGES  Plot the edges graphic.
%   SHOWEDGES(edges,vertex,varargin)  plot the edges graphic, in which the
%   edges is stored in E of n-by-1 matrix or cell, and vertex is
%   xy-coordinates of points,with plot properties pointed by varargin.
%
%   Example :
%       % Load a 2D triangulation and use the TriRep
%       % to construct a set of edges.
%       load trimesh2d;
%       % This loads triangulation tri and vertex coordinates x, y
%       trep = TriRep(tri, x,y)
%       e = edges(trep);
%       SHOWEDGES(E,[x,y],'r');
%
%       % Construct a triangulation using delaunay function
%       x = rand(30,1); y = rand(30,1);
%       TRI = delaunay(x,y);
%       e = edges(TRI);
%       SHOWEDGES(E,[x,y],'r');
%
%   See also edges, delaunay, TriRep, DelaunayTri.

%   Copyright 2012-2013 Jianbo Tang, CSU, China.
%   $Revision: 1.0.0.0 $
%

% checking inputs
if nargin<2||isempty(edges)||isempty(vertex)
    error('Not enough input arguments.');
end
if nargin<3||isempty(varargin{1})
    if iscell(edges)
        varargin = [];
    else
        varargin = {'m'};
    end
end

X = vertex(:,1);
Y = vertex(:,2);
Z = zeros(size(X));
option = true;
if size(vertex,2)>2
    Z = vertex(:,3);
    option = false;
end
clear vertex;

% plot figure
axes_nextplot = get(gca,'NextPlot');
if iscell(edges)
    nums = length(edges);  % subgraphs
    if isempty(varargin)
        c = hsv(nums);     % colors
        for i=1:nums
            E = edges{i};
            d = E(:,1:2)';
            h = plot3(X(d), Y(d), Z(d), 'Color', c(i,:));
            hold on;
        end
    else
        for i=1:nums
            E = edges{i};
            d = E(:,1:2)';
            h = plot3(X(d), Y(d), Z(d), varargin{:});
            hold on;
        end
    end
else
    d = edges(:,1:2)';
    h = plot3(X(d), Y(d), Z(d), varargin{:});
    hold on;
end
plot(X,Y,'ko','MarkerFaceColor','k','MarkerSize',3);
box on;
if option
    view(2);
end
set(gca, 'NextPlot', axes_nextplot);
if nargout
    varargout{1} = h;
end
end % showedges()


function L = edgesLength(E,X,Y)
%EDGESLENGTH  returns the length of edges in the triangulation.
%   L = EDGESLENGTH(E,x,y)  returns the length of edges in the
%   triangulation in a n-by-1 matrix, which n is edges numbers in the
%   triangulation, L(i) refers the length of edges E(i) returned by edges
%   function. x and y are coordinates used to construct the triangulation,
%   if L = EDGESLENGTH(E,X) called, X means the points x-y coordinate.
%
%   Example :
%       % Load a 2D triangulation and use the TriRep
%       % to construct a set of edges.
%       load trimesh2d;
%       % This loads triangulation tri and vertex coordinates x, y
%       trep = TriRep(tri, x,y)
%       e = edges(trep);
%       L = edgesLength(e,x,y);
%
%       % Construct a triangulation using delaunay function
%       x = rand(30,1); y = rand(30,1);
%       TRI = delaunay(x,y);
%       e = edges(TRI);
%       L = edgesLength(e,[x,y]);
%
%   See also edges, delaunay, TriRep, DelaunayTri.

%   Copyright 2012-2013 Jianbo Tang, CSU, China.
%   $Revision: 1.0.0.0 $
%

if nargin<2||isempty(E)||isempty(X)
    error('Not enough input arguments.');
end
if nargin==2
    % X<x,y>
    if size(X,2)<2
        error('Input argumets have invalid dimension.');
    end
    DX = X(E(:,1),1) - X(E(:,2),1);
    DY = X(E(:,1),2) - X(E(:,2),2);
    L  = sqrt(DX.^2+DY.^2);
elseif nargin==3
    if isempty(Y)
        error('Input argumets have invalid dimension.');
    end
    X = X(:);
    Y = Y(:);
    if length(X)~=length(Y)
        error('The latest inputs are not same dimension.');
    end
    DX = X(E(:,1)) - X(E(:,2));
    DY = Y(E(:,1)) - Y(E(:,2));
    L  = sqrt(DX.^2+DY.^2);
end % if
end % function: edgesLength()


