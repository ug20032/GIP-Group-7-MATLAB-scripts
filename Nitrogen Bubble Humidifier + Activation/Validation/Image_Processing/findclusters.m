function [loc,indices] = findclusters(Im,threshold,h_threshold)
    
    Im_logical = Im >= threshold;
       
    % Find connected components (i.e. clusters)
    labelled_clusters = bwlabel(Im_logical);
    
    % Number of clusters
    n_clusters = max(labelled_clusters(:));
      
    loc = [NaN NaN];
    indices = cell(1, n_clusters);

    for i = 1:n_clusters
        [rows, cols] = find(labelled_clusters == i);
        % Store indices in the cell array
        indices{i} = [rows, cols];

        ycentre = mean(rows);
        xcentre = mean(cols);
        loc(i, :) = [xcentre, ycentre]; % centre of clusters (indexed)

        if ycentre < h_threshold % remove 'bubble' if above surface
            loc(i, :) = [NaN NaN];
            indices{i} = [NaN];
        end
    end
end