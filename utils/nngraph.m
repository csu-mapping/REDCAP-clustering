function nngraph(geodata,connect,islabel)
%NNGRAPH  Show adjacent neighbor connect graph.
%   Syntax:
%      nngraph(geodata,connect,islabel)
%      nngraph(geodata,knn,islabel)
%   Input£º
%      geodata -- geographical coordinate <x,y>£»
%      connect -- spatial contiguity relatonshp, <cell npts-by-1>£»
%      islabel -- if show the node index on the graph£»
%   Output:
%      figure;
%
%   See also PLOT
%
%   Copyright 2012 Tang Jianbo, China.
%   This code may be freely used and distributed, so long as it maintains
%   this copyright line.
%   $Revision: 1.0 $     $Date: 2013/01/01 14:12:20 $

if nargin<2||isempty(geodata)
    error('  Not enough input arguments.');
end
if nargin<3||isempty(islabel)
    islabel = 'yes';
end

marker         = 'o';
markersize     = 12;
markercolor    = [0.5,0.5,0.5];
linewidth      = 2;
linecolor      = 'm';
if size(geodata,1)>20
    marker     = 'k.';
    markersize = 6;
    linewidth  = 1;
    islabel    = 'no';
end

if ~iscell(connect)
    connect = num2cell(connect, 2); % knn matrix
end
% show figures
figure;
hold on;
npts = size(geodata,1);
x = geodata(:,1);
y = geodata(:,2);
switch lower(islabel)
    case {'no', 0}
        for i=1:npts
            d = [i*ones(length(connect{i})); connect{i}];
            plot(x(d),y(d),'Color',linecolor,'LineWidth',linewidth);
        end
        plot(x,y,marker,'MarkerSize',markersize,'MarkerFaceColor',markercolor);
    otherwise
        for i=1:npts
            d = [i*ones(length(connect{i})); connect{i}];
            plot(x(d),y(d),'Color',linecolor,'LineWidth',linewidth);
        end
        plot(x,y,marker,'MarkerSize',markersize,'MarkerFaceColor',markercolor);
        for i=1:npts
            text(x(i),y(i),num2str(i),'HorizontalAlignment','center');
        end
end
box on;
end % function

