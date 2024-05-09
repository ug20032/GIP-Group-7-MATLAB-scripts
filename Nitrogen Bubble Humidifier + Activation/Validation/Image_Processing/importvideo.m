function [images_gray,images_rgb,video] = importvideo(filePath)

    disp('--- Importing video ---')
    video = VideoReader(filePath);
    images_gray = uint8(zeros(1280,720,video.NumFrames));
    images_rgb = uint8(zeros(1280,720,3,700));
    
    i = 1;
    tic
    disp('--- Converting video to images ---')
    while hasFrame(video)
    
        frame = readFrame(video);
        images_gray(:,:,i) = im2gray(frame); % Convert to grayscale
        images_rgb(:,:,:,i) = frame;
        
        i = i + 1;
        clockTime = toc;
        if mod(i,50) == 0
            disp([num2str(i/video.NumFrames*100)  '% complete, ' num2str(clockTime) 's elapsed'])
        end
    end
end