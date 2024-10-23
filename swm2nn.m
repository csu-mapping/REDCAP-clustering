function connect = swm2nn(spatial_weight_matrix)
% spatial_weight_matrix 
%     1 2   - [souceIndex, spatialLinkedIndex]
%     1 4
%     2 1
%     ...
%

n = max(max(spatial_weight_matrix));
connect = cell(n,1);

for i=1:n
    connect{i} = spatial_weight_matrix(spatial_weight_matrix(:,1)==i, 2)';
end

end 


