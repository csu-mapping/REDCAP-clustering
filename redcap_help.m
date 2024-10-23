%function [RegionIDX,Z,LinkageIDX,connect] = redcap(attrdata,connect,k,method,strategy)
%REDCAP  regionalization with dynamically constrained agglomerative
%        clustering and partitioning.
%  Syntax:
%     [RegionIDX,Z,LinkageIDX] = REDCAP(data,connect,k,method,strategy)
%       ...
%  Input��
%     attrdata - thematic dataset, <nxp double> includes n objects which with
%     p variables each;
%     connect - spatial contiguity cell matrix,<nx1 cell>, row j refers
%     the neighbors of point j in data;
%     k - max number of regions in partioning part;
%     method - agglomerative clustering mehod,e.g.,
%              Method       Description
%              'single'	    Shortest distance
%              'average'	Unweighted average distance (UPGMA)
%              'complete'	Furthest distance
%              'centroid'	Centroid distance (UPGMC), appropriate for
%              Euclidean distances only.
%              'ward'	    Inner squared distance (minimum variance
%              algorithm), appropriate for Euclidean distances only.
%     strategy - two constraining strategies,e.g.,'first','full'.
%  Output��
%     RegionIDX - cluster id of objects coincident the level<Nxregion_num>
%     Z - spatial congiuous tree;
%     LinkageIDX - cluster id of points in first agglomerative part.
%  Reference:
%     D. Guo (2008): Regionalization with dynamically constrained
%     agglomerative clustering and partitioning (REDCAP),
%     International Journal of Geographical Information
%     Science,22:7,801-823.
%
%  See also STRUCT2CELL, ISMEMBER, UNIQUE
%
%  Copyright 2013, Tang Jianbo, CSU, China.
%  This code may be freely used and distributed, so long as it maintains
%  this copyright line.
%  $Revision: 1.0 $     $Date: 2013/03/16 22:30:20 $

% This algorithm consists two parts,agglomerative clustering part and
% partioning part.
% This function includes six methods,coming from the combination of three
% linkage based clustering methods(SLK,ALK,CLK,WLK) and two constrained
% strategies(first-order,full-order).


