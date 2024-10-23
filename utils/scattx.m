function scattx(data,IDX,marker,markersize,colors,centerflag)
%SCATTX   display point clusters.
%   Synatx:
%         SCATTX(data,IDX,marker,markersize,centerflag)
%   Input:
%         data       - spatial dataset in which each row refers a point
%         IDX        - cluster id of every event in data,{[]}
%         marker     - makers to reprensent clusters,{12 types}
%         markersize - marker size,{6}
%         colors     - marker color,{16 types}
%         centerflag - show  the clusters center or not,{'off'}
%
%   See also  SCATTER, MAPSHOW, GEOSHOW, PLOTM
%
%   Copyright 2012 Tang Jianbo, CSU.
%   This code may be freely used and distributed, so long as it maintains
%   this copyright line.
%   Version 1.0, Date 2012-10-20 20:10:00.
%

%% Input arguments checking
if nargin<1||isempty(data)       % check - data
    error('Not enough input arguments.');
end
[npts,dims] = size(data);
if dims<2
    disp('Warning:');
    disp('   The first input argument should be more than 2 columns,');
    disp('   Each row refers to one spatial point in the dataset.');
    disp('   If you want to show one dimensional data to use ''plot(x)''.');
    return;
end

if nargin<2||isempty(IDX)        % check - IDX
    IDX = ones(npts,1);
else
    IDX = IDX(:);
    if npts ~= length(IDX)
        error(['"',inputname(1),'" and "',inputname(2),'" are not same length.']);
    end
end

if nargin<3||isempty(marker)    % check - marker
    marker = {'.','o','s','^','*','v','d','<','p','>','h','+'};  % 12 marker type
end
noise_marker = {'x'};

if nargin<4||isempty(markersize) % check - markersize
    markersize = 5;
end

if nargin<5||isempty(colors)
    % define color matrix
    colors = {   
                   [1,      0.502,  0.502];
                   [0,      0,      0];            % black
                   [0,      0,      1];            % blue
                   [1,  	0,  	0];            % red  
                   [0,      1,  	0];            % green
                   [1,      1,      0];            % yellow
                   [0,      1,      1];            % cyan
                   [1,      0,      1];            % magenta
                   [0.502,  0,      0];
                   [0,  	0.502,	0];            % darkgreen
                   [0.502,  0,      1];
                   [0.502,  0,  0.502];
                   [0,      1,      0.502];
                   [0.502,  0.502,  0]; 
                   [0,      0.502,  0.502];
                   [1,      0.502,  0];
              };                 % 16 color type
end
noise_color = [0.7, 0.7, 0.7]; % gray


if nargin<6||isempty(centerflag) % check - centerflag
    centerflag = 'off';  % do not show cluster centers on the figure
end

% cluster
EventsIDs = (1:npts)';
noise   = EventsIDs(IDX==0); % noise points indexs
clusterid = unique(IDX(IDX~=0)); % cluster id
cluster_nums = length(clusterid);

BL = (cluster_nums>(12*16));
if(BL)
    colors = num2cell(hsv(ceil(cluster_nums/12)),2);
    % when cluster number larger than 12*16, we caculate the color needed.
end
markernums = length(marker);
colornums  = length(colors);

%% plot figure
if isequal(get(gca,'NextPlot'),'replace')
    % find if the current figure is hold off;
    hold_on = false;
else
    hold_on = true;
end

switch lower(centerflag)
    case {'off','no','n',0}
        % doesn't show the cluster centers
        if(dims==2)
            % 2-D scatter plot
            if ~isempty(noise)
                plot(data(noise,1),data(noise,2),noise_marker{1},'MarkerSize',markersize,'MarkerFaceColor',noise_color,'color',noise_color);
                hold on;
            end
            for i = 1:cluster_nums
                X = EventsIDs(IDX==clusterid(i));
                % marker_id = rem(i,markernums)+1;
                % color_id  = ceil(i/markernums);
                color_id  = rem(i,colornums)+1;
                marker_id = rem(ceil(i/colornums),markernums)+1;
                plot(data(X,1),data(X,2),marker{marker_id},'MarkerSize',markersize,'MarkerFaceColor',colors{color_id},'color',colors{color_id});
                hold on;
            end
        elseif(dims>=3)
            % 3-D scatter plot
            if ~isempty(noise)
                plot3(data(noise,1),data(noise,2),data(noise,3),noise_marker{1},'MarkerSize',markersize,'MarkerFaceColor',noise_color,'color',noise_color);
                hold on;
            end
            for i=1:cluster_nums
                X = EventsIDs(IDX==clusterid(i));
                % marker_id = rem(i,markernums)+1;
                % color_id  = ceil(i/markernums);
                color_id  = rem(i,colornums)+1;
                marker_id = rem(ceil(i/colornums),markernums)+1;
                plot3(data(X,1),data(X,2),data(X,3),marker{marker_id},'MarkerSize',markersize,'MarkerFaceColor',colors{color_id},'color',colors{color_id});
                hold on;
            end
            view(-32,30);
            % zlabel('Time');
        else
            error('  The first paramter is more than 3 columns.');
        end % if
        
    otherwise       % case {'on','yes','y',1}
        % show cluster centers
        if(dims==2)
            % 2-D scatter plot
            if ~isempty(noise)
                plot(data(noise,1),data(noise,2),noise_marker{1},'MarkerSize',markersize,'MarkerFaceColor',noise_color,'color',noise_color);
                hold on;
            end
            for i = 1:cluster_nums
                X = EventsIDs(IDX==clusterid(i));
                % marker_id = rem(i,markernums)+1;
                % color_id  = ceil(i/markernums);
                color_id  = rem(i,colornums)+1;
                marker_id = rem(ceil(i/colornums),markernums)+1;
                plot(data(X,1),data(X,2),marker{marker_id},'MarkerSize',markersize,'MarkerFaceColor',colors{color_id},'color',colors{color_id});
                hold on;
                % plot cluster center
                center_x = mean(data(X,1));
                center_y = mean(data(X,2));
                text(center_x,center_y,num2str(i));
            end
        elseif(dims>=3)
            % 3-D scatter plot
            if ~isempty(noise)
                plot3(data(noise,1),data(noise,2),data(noise,3),noise_marker{1},'MarkerSize',markersize,'MarkerFaceColor',noise_color,'color',noise_color);
                hold on;
            end
            for i=1:cluster_nums
                X = EventsIDs(IDX==clusterid(i));
                % marker_id = rem(i,markernums)+1;
                % color_id  = ceil(i/markernums);
                color_id  = rem(i,colornums)+1;
                marker_id = rem(ceil(i/colornums),markernums)+1;
                plot3(data(X,1),data(X,2),data(X,3),marker{marker_id},'MarkerSize',markersize,'MarkerFaceColor',colors{color_id},'color',colors{color_id});
                hold on;
                % plot cluster center
                center_x = mean(data(X,1));
                center_y = mean(data(X,2));
                center_z = mean(data(X,3));
                text(center_x,center_y,center_z,num2str(i));
            end
            view(-32,30);
            % zlabel('Time');
        else
            error('  The first paramter is more than 3 columns.');
        end % if
end % switch()
%
% set(gcf,'color','w');
% xlabel('X');
% ylabel('Y');
box on;
% recover the original figure hold status
if hold_on
    hold on;
else
    hold off;
end

end % // scattx()

