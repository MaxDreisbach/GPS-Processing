clear all
clc
close all

Int_corr = 4;
Int_corr_BG = 4;
PLOT = false;
PLOT2 = false;
PLOT = true;
PLOT2 = true;
RUN_ON_CMD = true;
offset = 1; %offset of image from ground (1px necessary to get closed mask)
ground = 1;
rot_angle = 0.65; %tilt angle of substrate in degrees

%% Load Correction-Matrix
C_dat = load('RGB_Correction_Matrix_measured.mat');
C = cell2mat(struct2cell(C_dat));
%
%% Source: image path, image name
source_dir = '../Data/15032023_water_PDMS_0 degree-1_C001H001S0001/';
out = regexp(source_dir,'/','split');
out_dir = ['../Preprocessed_dil/',cell2mat(out(end-1)),'_preprocessed/'];

%% Get first image 
images = dir([source_dir,'*.tif']);  % Vector of all image filenames
Imax = imread([source_dir,images(1).name]);
Imax = double(Imax)/2^12;
[Imax] = Spectral_reconstruction(Imax,C,Int_corr_BG);
IBmax = Imax(:,:,3);
IBmax = imrotate(IBmax,rot_angle,'crop');
IBmax = IBmax(40:end-20,:);


%% Background
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

for i = 1:length(images)
    if PLOT
        i=50;
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
    
    %% make the image square   
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
    
    %% Create mask
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
    
    %Morphological opening
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

%close matlab in command line mode
if RUN_ON_CMD
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


function [I_rec] = Spectral_reconstruction(I,C,Int_corr)
    %This method performs a spectral correction of the input image I with
    %the correction Matrix C
    
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

    %% Reconstruct RGB-image
    I_rec(:,:,1) = RED;
    I_rec(:,:,2) = GREEN;
    I_rec(:,:,3) = BLUE;
end