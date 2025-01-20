%% Performs spectral corrrection of glare-point shadowgraphy images and creates videos for the bottom and side views

clear all
clc
close all

Int_corr = 4;

%% Source: image path, image name
source_dir = '../';
out_dir = '../videos/';

if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end

%% Load Correction-Matrix
C_dat = load('RGB_Correction_Matrix_measured.mat');
C = cell2mat(struct2cell(C_dat))


%Iterate over subfolders
D = dir(source_dir);
for k = 3:length(D) 
    sfol = D(k).name
    fullpath = strcat(source_dir,sfol,'/')
    
    images = dir([fullpath,'*.tif']);  % Vector of all image filenames
    clearvars I I_rec red_ch green_ch blue_ch RED GREEN BLUE
    
    %% Initialize output video
    videoRGB = VideoWriter([out_dir strcat(sfol,'_RGB.avi')],'Motion JPEG AVI'); %create the video object
    videoR = VideoWriter([out_dir strcat(sfol,'_below.avi')],'Motion JPEG AVI'); 
    videoG = VideoWriter([out_dir strcat(sfol,'_GP.avi')],'Motion JPEG AVI'); 
    videoB = VideoWriter([out_dir strcat(sfol,'_side.avi')],'Motion JPEG AVI'); 
    v.Quality = 100;
    v.FrameRate = 30;
    open(videoRGB); %open the file for writing
    open(videoR); %open the file for writing
    open(videoG); %open the file for writing
    open(videoB); %open the file for writing
    

    for i = 1:length(images)
        imgName = images(i).name;   
        
        %% Split channels of RGB-image
        I = imread([fullpath,imgName]);
        I = double(I)/2^12;
        red_ch = I(:,:,1)*Int_corr;
        green_ch = I(:,:,2)*Int_corr;
        blue_ch = I(:,:,3)*Int_corr;

        %% Correct grayscale images of the color channels
        RED = (C(1,1) * red_ch + C(1,2) * green_ch + C(1,3) * blue_ch);
        GREEN = C(2,1) * red_ch + C(2,2) * green_ch + C(2,3) * blue_ch;
        BLUE = C(3,1) * red_ch + C(3,2) * green_ch + C(3,3) * blue_ch;

        %% Only keep positive values   
        RED = max(RED,0);
        GREEN = max(GREEN,0);
        BLUE = max(BLUE,0);
        
        %plotCorrection("red", RED, GREEN, BLUE, red_ch, green_ch, blue_ch)
        %% Scaling image intensity
        f_red = max(max(RED))/max(max(red_ch));
        f_green = max(max(GREEN))/max(max(green_ch));
        f_blue = max(max(BLUE))/max(max(blue_ch));
        RED = RED/f_red;
        GREEN = GREEN/f_green;
        BLUE = BLUE/f_blue;
        
        %% Clip values   
        RED = min(RED,1);
        GREEN = min(GREEN,1);
        BLUE = min(BLUE,1);
        
        %% Reconstruct RGB-image
        I_rec(:,:,1) = RED;
        I_rec(:,:,2) = GREEN;
        I_rec(:,:,3) = BLUE;
        
        %% Write video frames
        writeVideo(videoRGB,I_rec);
        writeVideo(videoR,RED);
        writeVideo(videoG,GREEN);
        writeVideo(videoB,BLUE);
        
       i
    end
    close(videoRGB); %close the file
    close(videoR); %close the file
    close(videoG); %close the file
    close(videoB); %close the file
end



function [] = plotCorrection(color, RED, GREEN, BLUE, red_ch, green_ch, blue_ch)

    if color == "red"
        figure('Name','Measured Data','NumberTitle','off');
        imshow(RED)
        figure('Name','Measured Data','NumberTitle','off');
        imshow(red_ch)
        figure('Name','Measured Data','NumberTitle','off');
        imshow((red_ch-RED))
    end
    
    if color == "green"
        figure('Name','Measured Data','NumberTitle','off');
        imshow(GREEN)
        figure('Name','Measured Data','NumberTitle','off');
        imshow(green_ch)
        figure('Name','Measured Data','NumberTitle','off');
        imshow((green_ch-GREEN))
    end
    
    if color == "blue"
        figure('Name','Measured Data','NumberTitle','off');
        imshow(BLUE)
        figure('Name','Measured Data','NumberTitle','off');
        imshow(blue_ch)
        figure('Name','Measured Data','NumberTitle','off');
        imshow((blue_ch-BLUE))
    end
end
