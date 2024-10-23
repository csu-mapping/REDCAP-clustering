### REDCAP clustering

This project contains the implementation of the spatial clustering method in paper "Regionalization with dynamically constrained agglomerative clustering and partitioning (REDCAP)" (https://www.tandfonline.com/doi/full/10.1080/13658810701674970). 

#### Usage:

1. Abstract your data, such as polyline, or polygon features into points. Simply, you can use their center points to represent them. Making a point data as follows:

```octave
points = [X Y Z1 Z2 Z3 ...]
```

where points is a (element numbers) row and (2+attribute numbers) columns- matrix. 
X,Y are the x- and y- coordinates of the points or sites, and Z1, Z2, ... are the survey variable values at the sites. Each row represent a point in space. The line numbers are the total number of points (or elements) in data.

2. Define the spatial neighborhoods of each point. You can use the 'spatialneighbors.m' function (in folder 'private') to generate the required spatial neighbors information as a variable named 'connect'.

3. Using the REDCAP to clustering the input data into different clusters. The usage of REDCAP can be found in 'redcap_help.m'.

  ```octave
  [RegionIDX,Z,LinkageIDX,connect] = redcap(attrdata,connect,k,method,strategy)
  ```

In general, the parameters of the clustering algorithm are set to as:

```octave
connect = spatialneighbors(points(:,1), points(:,2), 'delaunay', 1);
[RegionIDX,Z,LinkageIDX,connect] = redcap(points(:, 3:end), connect, size(points,1)-1, 'alk', 'first');
```

where ‘RegionIDX’ contains the hierarchical clustering results with size(points,1)-1 to 2 clusters. You can select the optimal clustering results with SIL index or you need.

4. To show the clustering result.

  ```octave
  IDX = RegionIDX(:, k);
  scattx(points(:,1:2), IDX); 
  ```

  