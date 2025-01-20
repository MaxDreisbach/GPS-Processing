clear all
clc
close all

%% Spectral Correction and Video Creation for Shadowgraphy Images
% This script performs spectral correction on glare-point shadowgraphy
% images and generates videos for different views (bottom, GP, and side).
%
% Dependencies:
% - Requires file with correction matrix

%% Configuration Parameters
Int_corr = 4; % Intensity correction factor for reconstructed images

%% Source and Output Directories
source_dir = '../';
out_dir = '../videos/';

if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end

%% Load Correction-Matrix
C_dat = load('RGB_Correction_Matrix_measured.mat');
C = cell2mat(struct2cell(C_dat))


%% Iterate Through Subfolders and Process Images
D = dir(source_dir);
for k = 3:length(D) 
    sfol = D(k).name
    fullpath = strcat(source_dir,sfol,'/')
    
    images = dir([fullpath,'*.tif']);  % Vector of all image filenames
    clearvars I I_rec red_ch green_ch blue_ch RED GREEN BLUE
    
    %% Initialize output videos
    videoRGB = VideoWriter([out_dir strcat(sfol,'_RGB.avi')],'Motion JPEG AVI'); %create the video object
    videoR = VideoWriter([out_dir strcat(sfol,'_below.avi')],'Motion JPEG AVI'); 
    videoG = VideoWriter([out_dir strcat(sfol,'_GP.avi')],'Motion JPEG AVI'); 
    videoB = VideoWriter([out_dir strcat(sfol,'_side.avi')],'Motion JPEG AVI'); 
    %v.Quality = 100;
    %v.FrameRate = 30;
    open(videoRGB); %open the file for writing
    open(videoR); %open the file for writing
    open(videoG); %open the file for writing
    open(videoB); %open the file for writing
    
    %% Process Each Image
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
     %close the video files
    close(videoRGB);
    close(videoR);
    close(videoG);
    close(videoB);
end

%% Function: plotCorrection
function [] = plotCorrection(color, RED, GREEN, BLUE, red_ch, green_ch, blue_ch)
% This helper function visualizes the results and errors of the 
% color correction for a specified color channel (red, green, or blue).
%
% Inputs:
%   color    - (string) Specifies the color channel to visualize. Valid
%              values are "red", "green", or "blue".
%   RED      - (2D matrix) Corrected intensity values for the red channel.
%   GREEN    - (2D matrix) Corrected intensity values for the green channel.
%   BLUE     - (2D matrix) Corrected intensity values for the blue channel.
%   red_ch   - (2D matrix) Original intensity values for the red channel.
%   green_ch - (2D matrix) Original intensity values for the green channel.
%   blue_ch  - (2D matrix) Original intensity values for the blue channel.

    if color == "red"
        figure('Name','input','NumberTitle','off');
        imshow(RED)
        figure('Name','color-corrected','NumberTitle','off');
        imshow(red_ch)
        figure('Name','error','NumberTitle','off');
        imshow((red_ch-RED))
    end
    
    if color == "green"
        figure('Name','input','NumberTitle','off');
        imshow(GREEN)
        figure('Name','color-corrected','NumberTitle','off');
        imshow(green_ch)
        figure('Name','error','NumberTitle','off');
        imshow((green_ch-GREEN))
    end
    
    if color == "blue"
        figure('Name','input','NumberTitle','off');
        imshow(BLUE)
        figure('Name','color-corrected','NumberTitle','off');
        imshow(blue_ch)
        figure('Name','error','NumberTitle','off');
        imshow((blue_ch-BLUE))
    end
end
