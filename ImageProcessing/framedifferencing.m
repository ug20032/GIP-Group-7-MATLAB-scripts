function [tracked_images, images_rgb] = framedifferencing(video, images_grey, framelist)

    disp('--- Tracking bubbles ---')
    tracked_images = uint8(zeros(1280,720));
    for i = 1:video.NumFrames - 1
        alpha = 0.4;
        delta = abs(images_grey(:,:,i+1)-images_grey(:,:,i)); % calculate change in pixels (i.e. movement)
        tracked_images(:,:,i+1) = alpha*tracked_images(:,:,i) + (1-alpha)*delta; % low pass filter
    
        if mod(i,1e2) == 0
            disp([num2str(i/video.NumFrames*100)  '% complete'])
        end
    
    end
    tracked_images = tracked_images;
   
    disp('--- Colouring images ---')
    % Shade the grey images with red to show frame differencing
    for i = framelist(1:end-1)
        images_rgb(:, :, 1, i) = images_grey(:,:,i) + 5*tracked_images(:,:,i+1);
        images_rgb(:, :, 2, i) = images_grey(:,:,i);
        images_rgb(:, :, 3, i) = images_grey(:,:,i);
    
        %disp([num2str(i) '/' num2str(max(framelist))])
    end

end