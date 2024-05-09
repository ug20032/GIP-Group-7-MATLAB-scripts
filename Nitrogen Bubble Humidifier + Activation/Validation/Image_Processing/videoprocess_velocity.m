clc
clear
close all

%% User selection
iter_num = '1b';

%% Select and scale video for specified iteration
filePath = ['Iteration' iter_num '_velocity.mp4'];

brightness_threshold = 40;
if iter_num == '1a'
    scale = 720 / 100; % 'scale' pixels = 1 mm
    height_threshold = 540; % pixel count at the surface of the water
elseif iter_num == '1b'
    scale = 720 / 100; % 'scale' pixels = 1 mm
    height_threshold = 510; % pixel count at the surface of the water
elseif iter_num == '2a'
    scale = 720 / 100; 
    height_threshold = 490; 
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
framelist = 1:video.NumFrames;
[tracked_images, images_rgb] = framedifferencing(video, images, framelist);

%% Locate the bubbles
formatting
figure
for i = framelist(1:end-1)

    [loc,indices] = findclusters(tracked_images(:,:,i),brightness_threshold,height_threshold);
    if isempty(loc)
        z_min(i) = NaN;
    else
        z_min(i) = min(loc(:,2));
        z_max(i) = max(loc(:,2));
    end

    % show the animation
    clf
    imshow(images_rgb(:,:,:,i))
    hold on
    plot(loc(:,1), loc(:,2), 'o', 'MarkerFaceColor', colour(1,:), 'MarkerEdgeColor', colour(1,:), 'MarkerSize', 2);  
    %yline(height_threshold,'w')
    drawnow

end

%% Calculate the bubble velocity

% Image processing to find mean velocity
% start the position at the datum (z = 0)
datum = max(z_max);
z_c = -(z_min - datum) /scale;  
z_c(isnan(z_c)) = 0;

fps = 480; % frames per second
t = (1:length(z_c))/fps;
[vel,~,~] = linearfit(t,z_c);

figure
hold on
scatter(t*1e3, z_c,'MarkerEdgeColor','k','Marker','.')
plot(t*1e3, vel*(t*1e3),'Color',colour(1,:))
xlabel('Time (ms)')
ylabel('Y-position (mm)')
legend('$\max(z_c)$',['Velocity: ' num2str(round(vel,3)) 'm/s'],'Interpreter','latex','Location','SE')
ylim([0 100])
disp(['Mean velocity: ' num2str(round(vel*1e2,1)) 'cm/s'])


%% Functions

function [vel,yfit,err] = linearfit(t,y_pos)

    % Linear regression fitting
    vel = 0;
    yfit = vel * (t*1e3);
    err(1) = inf;
    err(2) = sqrt(sum((yfit - y_pos).^2));
    i = 1;
    while err(i+1) < err(i)
        i = i + 1;
        vel = vel + 1e-4;
        yfit = vel * (t*1e3);
        err(i+1) = sqrt(sum((yfit - y_pos).^2));
    end

end