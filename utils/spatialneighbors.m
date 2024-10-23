function [connect,nums] = spatialneighbors(x, y, option, k, varargin)
%spatialneighbors  Get connected neighbors of points.
%   Syntax:
%     connect = spatialneighbors(x,y)
%     connect = spatialneighbors(x,y,k)
%     connect = spatialneighbors(x,y,'delaunay',k)
%     connect = spatialneighbors(x,y,'cut_delaunay',a)
%     connect = spatialneighbors(x,y,'global',k)
%     connect = spatialneighbors(x,y,'local',k)
%     connect = spatialneighbors(x,y,'knn',k)
%     connect = spatialneighbors(x,y,'eps',eps);
%     connect = spatialneighbors(Eg,npts,k)
%   Input:
%     * x -- x coordinate;
%     * y -- y coordinate;
%     * k -- k order for 'delaunay' model or k nearest number for 'knn'
%       model,{1};
%     * option -- model, {'eps'}|'knn'|'delaunay'|'global'|'local'|'cut_delaunay',
%       (1) 'eps': neighbors in the epslion neighborhood of a point;
%       (2) 'knn': k nearest neighbors;
%       (3) 'delaunay': neighbors who are connected to a point in a delaunay
%          trianglation with path less than k;
%       (4) 'global': neighbors who are connected to a point in a global 
%          constranted delaunay trianglation with path less than k;
%       (5) 'local': neighbors who are connected to a point in a global and 
%          local constranted delaunay trianglation with path less than k;
%       (6) 'cut_delaunay': cut the edge whose length exceeds the
%       designated threshold(function of the mean).
%     * alpha -- constant, varargin{1};
%     * beta  -- constant, varargin{2}
%   Output:
%     connect -- cell array, connect{i} is neighbors of point(i), which
%     contains the index of neigbor points of point(i);
%     nums -- number of members in each point neighbor;
%   Example:
%{
       x = rand(50,1);
       y = rand(50,1);
       d_connect = neighbours(x,y,'delaunay',1);
       k_connect = neighbours(x,y,'knn',6);
%}
%   See also kOrderEdges, trigraph, globalcut, localcut, mst

%   Copyright 2012 Tang Jianbo, China.
%   This code may be freely used and distributed, so long as it maintains
%   this copyright line.
%   $Revision: 1.0 $     $Date: 2013/01/01 14:12:20 $

if nargin<2||isempty(x)||isempty(y)
    error('Not enough input arguments.');
end
if nargin<3||isempty(option)
    option = 'eps';
end
if nargin<4||isempty(k)||ischar(k)
    k = 1;
end

if size(x,2)>1 && isscalar(y)
    % Inputs:  edges, npts, k, [option]
    if nargin<3
        k = 1;
    else
        k = option;
    end
    Eg = x;
    Eg = Eg(:,1:2);
    npts = y;
    connect = cell(npts,1);
    for i=1:npts
        nn = kOrderEdges(Eg, i, k);
        nn = unique(nn(:));
        connect{i} = (nn(nn~=i))';
    end
else
    if size(x,2)>=2
        error('The first input should be a vector(x-coordinate).');
    end
    % Inputs:  x    , y   , k, [option]
    x = x(:);
    y = y(:);
    if length(x)~=length(y)
        error('The first two inputs are not the same length.');
    end
    npts = length(x);
    connect = cell(npts,1);
    switch lower(option)
        case {'eps','epslion'}
            points = [x,y];
            pointids = 1:npts;
            for i=1:npts
                dis = distx(points(i,:), points);
                ids = pointids(dis<k);
                connect{i} = ids(ids~=i);
                clear ids;
            end
        case {'delaunay', 'del', 'tri'}
            if nargin<5||isempty(varargin{1})
                alpha = [];
            else
                alpha = varargin{1};
            end
            if npts<3
                error('Too less point.');
            else
                dt = DelaunayTri(x,y);
                Eg = edges(dt);
                Eg = Eg(:,1:2);
                if ~isempty(alpha)
                    L  = sqrt((x(Eg(:,1))-x(Eg(:,2))).^2+(y(Eg(:,1))-y(Eg(:,2))).^2);
                    Eg = Eg(L<=(mean(L)+alpha*std(L)),:);
                    clear L;
                end
                clear dt;
                for i=1:npts
                    nn = kOrderEdges(Eg, i, k);
                    nn = unique(nn(:));
                    connect{i} = (nn(nn~=i))';
                end
            end
        case {'global', 'global-cut'}
            % alpha = 1;
            if nargin<5||isempty(varargin{1})
                alpha = 1;
            else
                alpha = varargin{1};
            end
            G = trigraph(x,y);
            Eg = cell2mat(globalcut(G, alpha, 0, 3));
            clear G;
            Eg = Eg(:,1:2);
            for i=1:npts
                nn = kOrderEdges(Eg, i, k);
                nn = unique(nn(:));
                connect{i} = (nn(nn~=i))';
            end
        case {'local',  'local-cut'}
            % alpha = 1;
            % beta  = 1.5;
            if nargin<5||isempty(varargin{1})
                alpha = 1;
            else
                alpha = varargin{1};
            end
            if nargin<6||isempty(varargin{2})
                beta = 1;
            else
                beta = varargin{2};
            end
            G = trigraph(x,y);
            Eg = globalcut(G, alpha, 0, 3);
            Eg = cell2mat(localcut(Eg, beta, 0, 3));
            Eg = Eg(:,1:2);
            for i=1:npts
                nn = kOrderEdges(Eg, i, k);
                nn = unique(nn(:));
                connect{i} = (nn(nn~=i))';
            end
        case {'knn', 'nearest'}
            if npts<k
                error('Too less point.');
            else
                knn = knn_search([x,y],[],k);
                connect = num2cell(knn, 2);
            end
        case {'mst','minimum spanning tree'}
            Eg = mst(x,y,'off');
            Eg = Eg(:,1:2);
            npts = length(y);
            connect = cell(npts,1);
            for i=1:npts
                nn = kOrderEdges(Eg, i, k);
                nn = unique(nn(:));
                connect{i} = (nn(nn~=i))';
            end
        case {'cut_del','del_cut','cut_delaunay'}
            Egs = delaunayx(x, y, k, 0);
            connect = spatialneighbors(Egs,length(x),1);
        otherwise
            disp('   Input option name is not found.');
            disp('   set option as bellow:');
            disp('   (1) eps: epslion neighbors of a point.');
            disp('   (2) knn: k nearest neighbors.');
            disp('   (3) delaunay: neighbors who are connected to');
            disp('       a point in a delaunay trianglation with');
            disp('       path less than k.');
            disp('   (4) global: neighbors who are connected to a point');
            disp('       in a global constranted delaunay trianglation');
            disp('       with path less than k.');
            disp('   (5) local: neighbors who are connected to a point');
            disp('       in a global and local constranted delaunay ');
            disp('       trianglation with path less than k.');
            disp('   (6) cut_delaunay: cut the edge whose length');
            disp('       exceeds the designated threshold(constant multiple');
            disp('       the mean of all edge lengths).');
    end % switch
end % if
if nargout==2
    nums = cellfun(@length, connect);
end
end % function



%% parseinputs
function [Q,R,K,fident] = parseinputs(varargin)
% Check input and output
error(nargchk(1,3,nargin));

Q=varargin{1};

if nargin<2
    R=Q;
    fident = true;
else
    fident = false;
    R=varargin{2};
end

if isempty(R)
    fident = true;
    R=Q;
end

if ~fident
    fident = isequal(Q,R);
end

if nargin<3
    K=1;
else
    K=varargin{3};
end
end % parseinputs()

%% knn_search
function [IDX,D] = knn_search(varargin)
% KNNSEARCH   Linear k-nearest neighbor (KNN) search
% IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by eahc row of Q (m x d array).
% The results are stored in the (m x K) index array, IDX.
% IDX = knnsearch(Q,R) takes the default value K=1.
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.
%
% See also, kdtree, nnsearch, delaunary, dsearch
%
% By Yi Cao at Cranfield University on 25 March 2008
%

% Check inputs
[Q,R,K,fident] = parseinputs(varargin{:});

% Check outputs
error(nargoutchk(0,2,nargout));

% C2 = sum(C.*C,2)';
[N,M] = size(Q);
L=size(R,1);
IDX = zeros(N,K);
D = IDX;

if K==1
    % Loop for each query point
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [D(k),IDX(k)]=min(d);
    end
else
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [s,t]=sort(d);
        IDX(k,:)=t(1:K);
        D(k,:)=s(1:K);
    end
end
if nargout>1
    D=sqrt(D);
end
end % knnsearch()

%% distx
function D = distx(point,pointset)
%DISTX   compute euclidean distance between one point and a point set.
if nargin<2||isempty(point)||isempty(pointset)
    error('Not enough input arguments.');
end
if size(point,2)~=size(pointset,2)
    error('The two arguments are not the same dimension.');
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


function E = kOrderEdges(graph_edges,vertex_id,k)
%KORDEREDGES    Searching the k order connected edges with the point VI.
%   E = kOrderEdges(GE,VI,K)  Searching the k order connected edges with
%   the point VI, GE is a n-by-2(3) matrix and GE(i) refers a edge in
%   delaunay graph, K is the order number.
%
%   Example :
%{
      % Load a 2D triangulation and use the TriRep
      % to construct a set of edges.
      load trimesh2d;
      % This loads triangulation tri and vertex coordinates x, y
      trep = TriRep(tri, x,y)
      e = edges(trep);
      knnE = kOrderEdges(e, 2, 1);

      % Construct a triangulation using delaunay function
      x = rand(30,1); y = rand(30,1);
      TRI = delaunay(x,y);
      e = edges(TRI);
      knnE = kOrderEdges(e, 2, 2);
%}
%   See also edges, delaunay, TriRep, DelaunayTri.

%   Copyright 2012-2013 Jianbo Tang, CSU, China.
%   $Revision: 1.0.0.0 $

% Searching the k-connected edges in the Delaunay network
if nargin<2||isempty(graph_edges)||isempty(vertex_id)
    error('Not enough input arguments.');
end
if nargin<3
    k = 1;
end

% When K=1
E = unique([graph_edges(graph_edges(:,1)==vertex_id,:); ...
    graph_edges(graph_edges(:,2)==vertex_id,:)], 'rows');

% When K>1
for i=1:(k-1)
    PointID = unique([E(:,1); E(:,2)]);
    PointNums = length(PointID);
    for j=1:PointNums   
        % searching the directly connected edges with pointid(j)
        E = [E; graph_edges(graph_edges(:,1)==PointID(j),:); ...
            graph_edges(graph_edges(:,2)==PointID(j),:)]; 
    end
    E = unique(E,'rows');
end % for
% E = unique(E,'rows');
end %  kOrderEdges


function [G,X] = trigraph(varargin)
%TRIGRAPH  Delaunay triangular graph.
%  Syntax:
%     G = trigraph(X);
%     G = trigraph(x,y);
%     G = trigraph(...,c);
%     more detail search DelaunayTri for help;
%  Inputs:
%     X, (x,y) -- coordinates;
%     c        -- constraint egdes;
%  Output:
%     G -- graph,<nedge-by-3 double>,each row with respact to one
%     edge(start node, end node) and its length.
%
%  See also DelaunayTri
%

dt = DelaunayTri(varargin{:});
e  = edges(dt);
dx = dt.X(e(:,1),1) - dt.X(e(:,2),1);
dy = dt.X(e(:,1),2) - dt.X(e(:,2),2);
edge_length = sqrt(dx.^2+dy.^2);
clear dx dy;
G = [e, edge_length];

if nargout>1
    X = dt.X;
end
end % trigraph


function globalnet = globalcut(G, alpha, mvd, minpts)
%GLOBALCUT  global constraint delaunay triangulations.
%  Syntax:
%     gnet = globalcut(G, alpha);
%  Inputs:
%     G -- <n-by-3 double>,graphic network,e_start-e_end-e_length;
%     alpha -- constant;
%     mvd -- (optional argument);
%     minpts -- (optional argument); 
%  Output:
%     globalnet -- <m-by-1 cell>,delaunay network sub-graph with global
%     constrant.
%

if nargin<1||isempty(G)
    error('Not enough input arguments.');
end
if nargin<2||isempty(alpha)
    alpha = 1;
end
if nargin<3||isempty(mvd)
    mvd = 0; % minimum edge length threshold 
end
if nargin<4||isempty(minpts)
    minpts = 3; % minimum members within a graph
end

% subgraph
G = subgraph(G);
% global long edge cutting
globalnet = [];
for i=1:length(G)
    globalnet = [globalnet; global_constraint_cut(G{i},alpha,mvd,minpts)];
end
end % globalcut


%%
function NG = global_constraint_cut(G, alpha, mvd, minpts)
if size(G,1)<3||size(G,2)~=3
    error('Input argument is not valid.');
end

% global constranit
pointids = unique([G(:,1); G(:,2)]);
pointnum = length(pointids);

global_mean = mean(G(:,3));
global_std  = std(G(:,3));
edge_nums = size(G,1);
isOtherEdge = true(edge_nums,1);
G = [G, (1:edge_nums)'];

for i=1:pointnum
    id = pointids(i);
    edges = [G(G(:,1)==id,:); G(G(:,2)==id,:)];
    edgesnum = size(edges,1);
    if edgesnum>=1
        local_mean = mean(edges(:,3));
        threshold = global_mean+alpha*(global_mean/local_mean)*global_std;
        for j=1:edgesnum
            if (edges(j,3)>=threshold) && (edges(j,3)>=mvd)
                isOtherEdge(edges(j,4)) = false;
            end
        end
    end
end
clear pointids edges;

NG = G(isOtherEdge,1:3);
clear G isOtherEdge;
NG = subgraph(NG, minpts);

end % global_constraint_cut


function localnet = localcut(G, beta, mvd, minpts)
%LOCALCUT  global constraint delaunay triangulations.
%  Syntax:
%     gnet = localcut(G, alpha);
%  Inputs:
%     G -- <n-by-3 double>,graphic network,e_start-e_end-e_length;
%     beta -- constant;
%     mvd -- (optional argument);
%     minpts -- (optional argument); 
%  Output:
%     localnet -- <m-by-1 cell>,delaunay network sub-graph with local
%     constrant.
%

if nargin<1||isempty(G)
    error('Not enough input arguments.');
end
if nargin<2||isempty(beta)
    beta = 1;
end
if nargin<3||isempty(mvd)
    mvd = 0; % minimum edge length threshold 
end
if nargin<4||isempty(minpts)
    minpts = 3; % minimum members within a graph
end

% subgraph
G = subgraph(G);
% local long edge cutting
localnet = [];
for i=1:length(G)
    localnet = [localnet; local_constraint_cut(G{i},beta,mvd,minpts)];
end
end % localcut


%%
function net = local_constraint_cut(G,beta,mvd,minpts)
if size(G,1)<3||size(G,2)<3
    disp('   Edges in the graph are less than three.');
    net = G;
    return;
end

% local constranit
pointids = unique([G(:,1); G(:,2)]);
pointnum = length(pointids);
edge_nums = size(G,1);
G = [G, (1:edge_nums)'];

local_std = zeros(pointnum,1);
% mean_variation = 0;
knn_edges = cell(pointnum,1);

% compute local std
for i=1:pointnum
    % ponit id
    id = pointids(i);
    % directly connect edges
    directly_connect_edges = [G(G(:,1)==id,:); G(G(:,2)==id,:)];
    local_std(i) = std(directly_connect_edges(:,3));
    % mean_variation = mean_variation+std(directly_connect_edges(:,3));
    % searching two-order connected edges of point id
    nn = unique([directly_connect_edges(:,1); directly_connect_edges(:,2)]);
    nn = nn(nn~=id);
    knn_edges{i} = directly_connect_edges;
    for j=1:length(nn)
        nnid = nn(j);
        nn_dc_edges = [G(G(:,1)==nnid,:); G(G(:,2)==nnid,:)];
        knn_edges{i} = [knn_edges{i}; nn_dc_edges];
    end
    knn_edges{i} = unique(knn_edges{i},'rows');
end
clear directly_connect_edges  nn  nn_dc_edges  id nnid;

% local constranit
% mean_variation = mean_variation/pointnum;
mean_variation = mean(local_std);
clear local_std;
isOtherEdge = true(edge_nums,1);

for i=1:pointnum
    korder_edges = knn_edges{i};
    edgesnum = size(korder_edges,1);
    if edgesnum>=1
        korder_mean = mean(korder_edges(:,3));
        threshold = korder_mean + beta*mean_variation;
        for j=1:edgesnum
            if (korder_edges(j,3)>=threshold)&&(korder_edges(j,3)>=mvd)
                isOtherEdge(korder_edges(j,4)) = false;
            end
        end
    end
end
clear knn_edges pointids korder_edges;

net = G(isOtherEdge,1:3);
clear G isOtherEdge;
net = subgraph(net, minpts);

end % local_constraint_cut


function Z = mst(varargin)  
%MST  Constructing a eulicdean minimum spanning tree of scatter points.
%   Syntax:
%      Z = mst(edges)
%      Z = mst(x, y)
%      Z = mst(x, y, 'yes')
%   Input:
%      edges -- edges,<p1,p2,length>;
%      x -- x coordinate of points;
%      y -- y coordinate of points;
%      fig -- whether to show the minimum spanning tree figure;
%   Output:
%      Z -- edges of the minimum spanning tree,<m-by-3> in the format of
%          [point1,point2,weight].
%
%   See also  biograph, minspantree
%
%   Copyright 2012 Tang Jianbo, China.
%   This code may be freely used and distributed, so long as it maintains
%   this copyright line.
%   $Revision: 1.0 $     $Date: 2013/01/01 14:12:20 $

if nargin==1
    % Z = mst(edges)
    E = varargin{1};
    nMST = grMinSpanTree(E);
    Z = E(nMST,:);
    Z = unique(sortrows(Z,[1,2]),'rows');
elseif nargin==2
    % Z = mst(x, y)
    x = varargin{1};
    y = varargin{2};
    x = x(:);
    y = y(:);
    if length(x) ~= length(y)
        error('The inputs are not the same length.');
    end
    E = trigraph(x,y);
    nMST = grMinSpanTree(E);
    Z = E(nMST,:);
    Z = unique(sortrows(Z,[1,2]),'rows');
elseif nargin==3
    % Z = mst(x, y, 'yes')
    x = varargin{1};
    y = varargin{2};
    fig = varargin{3};
    x = x(:);
    y = y(:);
    if length(x) ~= length(y)
        error('The inputs are not the same length.');
    end
    E = trigraph(x,y);
    nMST = grMinSpanTree(E);
    Z = E(nMST,:);
    Z = unique(sortrows(Z,[1,2]),'rows');
    % figure
    switch lower(fig)
        case {'on', 'yes', 'y', 1}
            % figure
            d = Z(:,1:2)';
            plot(x(d), y(d), 'r');
            hold on;
            plot(x,y,'ko','MarkerFaceColor','k','MarkerSize',3);
            hold off;
        otherwise
    end
else
    error('Not enough or too many input arguments.');
end % if

end % MST



%% grMinSpanTree()
function nMST = grMinSpanTree(E)
% Function nMST=grMinSpanTree(E) solve 
% the minimal spanning tree problem for a connected graph.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of edge;
%     m - number of edges.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   nMST(n-1,1) - the list of the numbers of edges included 
%     in the minimal (weighted) spanning tree in the including order.
% Uses the greedy algorithm.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1
    error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

% The data preparation
En=[(1:m)',E]; % we add the numbers
En(:,2:3) = sort(En(:,2:3), 2); % edges on increase order
ln=find(En(:,2)==En(:,3)); % the loops numbers
En=En(setdiff((1:size(En,1))',ln),:); % we delete the loops
[~,iw]=sort(En(:,4)); % sort by weight
Ens=En(iw,:); % sorted edges
clear E En ln iw;
% We build the minimal spanning tree by the greedy algorithm 
[edges_nums, ndim] = size(Ens);
Emst = zeros(n-1, ndim);
Emst(1,:) = Ens(1,:); % 1st edge include to minimal spanning tree
count_mst = 1;
% Ens = Ens(2:end,:); % rested edges
start_edg = 2;
while (count_mst<n-1) && (start_edg<edges_nums)
    count_mst = count_mst + 1;
    Emst(count_mst,:) = Ens(start_edg,:); % we add next edge to spanning tree
    % Ens=Ens(2:end,:); % rested edges
    start_edg = start_edg + 1;
    if any( (Emst(count_mst,2)==Emst(1:count_mst-1,2)) & ...
            (Emst(count_mst,3)==Emst(1:count_mst-1,3)))|| ...
            IsCycle(Emst(1:count_mst,2:3))  % the multiple edge or cycle
        % Emst=Emst(1:end-1,:); % we delete the last added edge
        Emst(count_mst,:) = zeros(1,ndim);
        count_mst = count_mst - 1;
    end
end
nMST = Emst(1:(n-1),1); % numbers of edges
end % function

%% IsCycle()
function ic = IsCycle(E)
% true, if graph E have cycle
n=max(max(E)); % number of vertexes
A=zeros(n);
A((E(:,1)-1)*n+E(:,2))=1;
A=A+A'; % the connectivity matrix
p=sum(A); % the vertexes power
ic=false;
while any(p<=1), % we delete all tails
  nc=find(p>1); % rested vertexes
  if isempty(nc),
    return
  end
  A=A(nc,nc); % new connectivity matrix
  p=sum(A); % new powers
end
ic=true;
end % function

%% grValidation()
function [m,n,newE] = grValidation(E)
% The validation of array E - auxiliary function for GrTheory Toolbox.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if ~isnumeric(E),
  error('The array E must be numeric!') 
end
if ~isreal(E),
  error('The array E must be real!') 
end
se=size(E); % size of array E
if length(se)~=2,
  error('The array E must be 2D!') 
end
if (se(2)<2),
  error('The array E must have 2 or 3 columns!'), 
end
if ~all(all(E(:,1:2)>0)),
  error('1st and 2nd columns of the array E must be positive!')
end
if ~all(all((E(:,1:2)==round(E(:,1:2))))),
  error('1st and 2nd columns of the array E must be integer!')
end
m=se(1);
if se(2)<3, % not set the weight
  E(:,3)=1; % all weights =1
end
newE=E(:,1:3);
n=max(max(newE(:,1:2))); % number of vertexes
end % function

