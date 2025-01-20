clear all
clc
close all

%% Preprocess images for volumetric reconstruction
% This script processes glare-point shadowgraphy images of droplets for volumetric reconstruction.
% It applies spectral correction, crops and rotates images, and generates binary masks.
% The output includes preprocessed RGB images and binary mask files.
% 
% Usage:
% - Set `source_dir` to the directory containing input TIFF images.
% - Configure parameters like `rot_angle` for substrate tilt correction.
% - Enable `PLOT` and `PLOT2` flags for visualization of intermediate steps.
% - Run the script in MATLAB; set `RUN_ON_CMD` to true for automated command-line usage.
%
% Dependencies:
% - Requires `RGB_Correction_Matrix_measured.mat` for spectral correction.
% - Uses `Spectral_reconstruction` function to correct RGB cross-talk.
% - Relies on image processing functions for masking and morphological operations.
%
% Outputs:
% - Preprocessed images and binary masks saved to `out_dir`.

%% Configuration Parameters
Int_corr = 4; % Intensity correction factor for reconstructed images
Int_corr_BG = 4; % Intensity correction for background images
PLOT = true; % Enable visualization of processing steps
PLOT2 = true; % Enable detailed visualization
RUN_ON_CMD = true; % Close MATLAB after execution (for command-line usage)
offset = 1; % Image offset from ground (ensures mask closure)
ground = 1; % Ground pixel offset
rot_angle = 0.65; % Tilt angle of substrate in degrees

%% Load Correction-Matrix
C_dat = load('RGB_Correction_Matrix_measured.mat');
C = cell2mat(struct2cell(C_dat));

%% Source: image path, image name
source_dir = '../Data/15032023_water_PDMS_0 degree-1_C001H001S0001/';
out = regexp(source_dir,'/','split');
out_dir = ['../Preprocessed_dil/',cell2mat(out(end-1)),'_preprocessed/'];


%% Process First Image to Initialize
images = dir([source_dir,'*.tif']);  % Vector of all image filenames
Imax = imread([source_dir,images(1).name]);
Imax = double(Imax)/2^12;
[Imax] = Spectral_reconstruction(Imax,C,Int_corr_BG);
IBmax = Imax(:,:,3);
IBmax = imrotate(IBmax,rot_angle,'crop');
IBmax = IBmax(40:end-20,:);


%% Background Image Processing
% Iterate through initial images to extract background estimate.
for i = 1:5
    imgName = images(i).name;   
    I = imread([source_dir,imgName]);
    I = double(I)/2^12;
    I = imrotate(I,rot_angle,'crop');
    I = I(40:end-20,:,:); 
    [I_rec] = Spectral_reconstruction(I,C,Int_corr_BG);
    IB = I_rec(:,:,3);
    
    for i=1:size(IB,1)
        for j=1:size(IB,2)
            if IB(i,j)>IBmax(i,j)
                IBmax(i,j)=IB(i,j);
            end
        end
    end
end

% Calculate average background intensities
BG = Imax(1:end-100,:,:);
RED_BG = mean(mean(BG(:,:,1)))/Int_corr_BG*Int_corr;
GREEN_BG = mean(mean(BG(:,:,2)))/Int_corr_BG*Int_corr;
BLUE_BG = mean(mean(BG(:,:,3)))/Int_corr_BG*Int_corr;

%eliminate droplet from background image
IBmax(1:256,:) = BLUE_BG;

if PLOT2
    figure('Name','background')
    imshow(IBmax)
end

%% Process and Save All Images
for i = 1:length(images)
    if PLOT
        i=50; % Debug: Process and visualize only one image
    end
    imgName = images(i).name;  
    
    %% Split channels of RGB-image
    I = imread([source_dir,imgName]);
    I = double(I)/2^12;
    I = imrotate(I,rot_angle,'crop');
    %I = I(40:end-20,25:end-25,:);
    I = I(40:end-20,:,:); 
    [I] = Spectral_reconstruction(I,C,Int_corr);
    red_ch = I(:,:,1);
    green_ch = I(:,:,2);
    blue_ch = I(:,:,3);
    
    %% Create Square Image with Background Padding
    [height, width] = size(blue_ch);
    margin = abs(width-height);
    
    square_R = zeros(width,width)*RED_BG;
    square_G = ones(width,width)*GREEN_BG;
    square_B = ones(width,width)*BLUE_BG;
    square_R(margin+1-offset:width-offset,:) = red_ch;
    square_G(margin+1-offset:width-offset,:) = green_ch;
    square_B(margin+1-offset:width-offset,:) = blue_ch;
    square_B(end-offset:end,:) = 0;
    
    I_square = cat(3,square_R,square_G,square_B);
    
    %% Create Binary Mask
    BG_mask = imcomplement(imbinarize(IBmax));
    BW = imbinarize(blue_ch*4);
    BW2 = imcomplement(BW);
    BW3 = imfill(BW2,'holes');
    BW4 = BW3 - BG_mask;
    BW5 = imfill(BW4,'holes');
    Mask_d = double(BW5);
    
    if PLOT2
        figure('Name','BW')
        imshow(BW)
        figure('Name','BW2')
        imshow(BW2)
        figure('Name','BW5')
        imshow(BW5)
    end
    
    % Morphological opening for noise reduction
    se = strel('disk',4,4);
    Mask_d = imerode(Mask_d, se);
    Mask_d = imdilate(Mask_d, se);

    if PLOT2
        figure('Name','Mask_d')
        imshow(Mask_d)
    end
    
    % Dilate Mask in accordance with training data
    se = strel('disk',10,6);
    Mask_dil = imdilate(Mask_d, se);
    Mask_dil(end-offset-ground:end,:) = 0;
    Mask_d = Mask_dil;

    % Subtract structured ground
    Mask_d = Mask_d - BG_mask;
    Mask_d = imfill(Mask_d,'holes');

    if PLOT2
        figure('Name','Dilated Mask')
        imshow(Mask_d)
    end

    square_B_bg = zeros(width,width);
    square_B_bg(margin+1-offset:width-offset,:) = Mask_d;
    square_B_bg(end-offset:end,:) = 0;
    Mask_d = double(square_B_bg);    
    
    if PLOT2
        figure('Name','Mask holes filled')
        imshow(Mask_d)
        figure('Name','Image')
        imshow(I_square)
    end

    %% Rescale and save images
    cut = 200;
    top = 400;
    bottom = 000;
    left = 430;
    right = -30;
    I_square = I_square(cut*2+top:end-bottom,cut+left:end-cut-right,:);
    Mask_d = Mask_d(cut*2+top:end-bottom,cut+left:end-cut-right,:);
    I = imresize(I_square,[512 512]);
    M = imresize(Mask_d,[512 512]);
    
    if ~exist(out_dir)
       mkdir(out_dir) 
    end
    
        
    imwrite(I,[out_dir, imgName(1:end-4),'.png'])
    imwrite(M,[out_dir, imgName(1:end-4),'_mask.png'])
    %}
    if PLOT
        maskedRgbImage = bsxfun(@times, I_square, cast((Mask_d), 'like', I_square));
        figure('Name','Applied Mask')
        imshow(maskedRgbImage)
        return
    end
    i 
end


if RUN_ON_CMD
    %close matlab in command line mode
    id  = feature('getpid');
    if ispc
      cmd = sprintf('Taskkill /PID %d /F',id);
    elseif (ismac || isunix)
      cmd = sprintf('kill -9 %d',id);
    else
      disp('unknown operating system');
    end
    system(cmd);
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



function [I_rec] = Spectral_reconstruction(I,C,Int_corr)
    % This function adjusts the color intensities of an input RGB image
    % using a correction matrix to account for cross-talk. It scales
    % the corrected image channels to ensure consistency in intensity levels.
    %
    % Inputs:
    %   I        - (HxWx3 double) Input RGB image, normalized to [0, 1].
    %   C        - (3x3 double) Correction matrix for spectral correction.
    %              Each element represents the contribution of one light 
    %              source to an image channel.
    %   Int_corr - (double) Intensity correction factor
    %
    % Outputs:
    %   I_rec    - (HxWx3 double) Spectrally corrected RGB image, with intensity
    %              values adjusted and scaled to maintain consistency.
    
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

    %% Reconstruct the corrected RGB image
    I_rec(:,:,1) = RED;
    I_rec(:,:,2) = GREEN;
    I_rec(:,:,3) = BLUE;
end