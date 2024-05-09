function [d_bubble,loc2,indices2] = sizeclusters(loc,indices,threshold)

    % threshold = the pixel length threshold for a cluster

    % d_bubble = the distance between the two furthest points of the array
    % of cluster points (indices)
    % loc2 = same as loc but with rows removed for clusters that do not
    % meet threshold size
    
    % find the size of each cluster
    d_bubble = zeros(1,length(indices));
    for j = 1:length(indices) 
        d_bubble(j) = max(max(pdist2(indices{j}, indices{j})));
    end
    % remove all clusters which are outside the pixel size threhsold
    % for diameters
    d_max_logical = d_bubble > threshold(1) & d_bubble < threshold(2);
    d_bubble = d_bubble .* d_max_logical;
    d_bubble = d_bubble(d_bubble ~= 0);
    d_bubble(isnan(d_bubble)) = []; % remove NaNs

    % for centre of bubble
    loc(~d_max_logical, :) = NaN; 
    loc2 = loc(~any(isnan(loc), 2), :);
    
    % for all points in bubble
    indices2 = indices(d_max_logical);

end