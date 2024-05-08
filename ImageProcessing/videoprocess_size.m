clc
clear
close all
formatting
% Choose which iteration to select
iter_num = '2b';
filePath = ['Iteration' iter_num '_size.mp4'];

% brightness threshold for finding clusters
brightness_threshold = 7; 

% pixel distance threshold for sizing clusters
pixel_threshold = [40 200];

if iter_num == '1a'
    scale = 720 / 100; % 'scale' pixels = 1 mm
    height_threshold = 560; % pixel count at the surface of the water
elseif iter_num == '1b'
    scale = 720 / 100;
    height_threshold = 530;
elseif iter_num == '2a'
    scale = 720 / 100; 
    height_threshold = 510; 
elseif iter_num == '2b'
    scale = 720 / 100; 
    height_threshold = 550; 
elseif iter_num == '3a' 
    scale = 1280 / 167; 
    height_threshold = 530; 
elseif iter_num == '3b'
    scale = 1280 / (160 - 24); 
    height_threshold = 650;
end

% Import video
[images,images_orig, video] = importvideo(filePath);

% Track the bubbles
framelist = 1:25;
[tracked_images, images_rgb] = framedifferencing(video, images, framelist);

%% Locate and size the bubbles
figure
%for i = framelist(2:end-1)
for i = 17

    % find and size clusters (i.e. bubbles)
    [loc,indices] = findclusters(tracked_images(:,:,i),brightness_threshold,height_threshold);
    [d_max,loc,indices2] = sizeclusters(loc,indices,pixel_threshold); 

    % calculate & output the major semi-axis (r_a) 
    d_max_mm = d_max / scale;
    d_mean(i) = mean(d_max_mm); % mean across all bubbles in a single snapshot
    d_time_average(i) = mean(d_mean(2:end)); % mean across all bubbles over all snapshots (should converge)
    disp(['r_a = ' num2str(d_time_average(i)/2) ' mm, frames analysed: ' num2str(i) '/' num2str(framelist(end-1))])
 
    % show the animation
    clf
    imshow(images_rgb(:,:,:,i))
    hold on
    plot(loc(:,1), loc(:,2), 'o', 'MarkerFaceColor', colour(1,:), 'MarkerEdgeColor', colour(1,:), 'MarkerSize', 4); 
    %yline(height_threshold,'w')
    co_ords = bubbleboundary(indices2);
    drawnow

    print('frame_black8', '-dpng', '-r600');

end

%% Functions

function co_ords = bubbleboundary(indices2)

    % Plot boundary of bubble (convex hull)
    for ii = 1:length(indices2)
        co_ords = indices2{ii};
        k = convhull(co_ords(:,1), co_ords(:,2));
        plot(co_ords(k,2), co_ords(k,1), 'k-','Linewidth',0.5);
    end

end

