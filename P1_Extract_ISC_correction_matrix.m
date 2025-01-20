clear all
clc
close all

%% In-situ Color Correction: Image Processing Script
% This script constructs a transfer matrix for in-situ RGB correction using
% images captured under red, green, and blue illuminations. It calculates a
% correction matrix and applies it to a test image for spectral
% reconstruction. Corrected and original images are visualized and saved.
%
% Usage:
% - Set the `rpath` variable to the folder containing the image sequence of
% the experiments conducted with only red illumination. 
% Accordingly set `gpath` and `bpath`.


%% Source: image path, image name
% recordings with only red illumination
rpath = '17062022_Makrolon_hydrophobic_100percent_3cm_135deg_red\';
rfiles = dir([rpath,'*.tif']);
% green illumination only
gpath = '17062022_Makrolon_hydrophobic_100percent_3cm_135deg_green\';
gfiles = dir([gpath,'*.tif']);
% blue illumination only
bpath = '17062022_Makrolon_hydrophobic_100percent_3cm_blue\';
bfiles = dir([bpath,'*.tif']);

%% Get mean intensities of single LED images
for i = 1:length(rfiles)
  
Red = imread([rpath, rfiles(i).name]);
[Rr,Rg,Rb] = getIntensities(Red); %Red LED - red channel, Red LED - green channel, Red LED - blue channel

Green = imread([gpath, gfiles(i).name]);
[Gr,Gg,Gb] = getIntensities(Green); %Green LED - red, green, blue channel

Blue = imread([bpath, bfiles(i).name]);
[Br,Bg,Bb] = getIntensities(Blue); %Blue LED - red, green, blue channel

%% Construct transfer matrix
T = [Rr Gr Br; Rg Gg Bg; Rb Gb Bb];
T_Stack(:,:,i) = T;

%Create composite image for visualisation
Comp = cat(3,Red(:,:,1)*60, Green(:,:,2)*60, Blue(:,:,3)*60);
figure
imshow(Comp);
end

%% Averaging and calculatio of correction matrix
T_mean = mean(T_Stack,3);
C = inv(T_mean);
Cn = C/(max(max(C)));

filename = 'RGB_Correction_Matrix_measured.mat';
save(filename,'Cn')


%% Testing - Spectral reconstruction of recordings
testpath = "12052022_Makrolon_hydrophobic_45_115_50percent_high_01\12052022_Makrolon_hydrophobic_45_115_50percent_high_01000029.tif";
I = imread(testpath);

C = Cn;

red_ch = double(I(:,:,1));
green_ch = double(I(:,:,2));
blue_ch = double(I(:,:,3));

RED = (C(1,1) * red_ch + C(1,2) * green_ch + C(1,3) * blue_ch);
GREEN = C(2,1) * red_ch + C(2,2) * green_ch + C(2,3) * blue_ch;
BLUE = C(3,1) * red_ch + C(3,2) * green_ch + C(3,3) * blue_ch;

figure('Name', 'RGB')
imshow(double(I)/max(max(blue_ch))*10)

figure('Name', 'Red corrected')
imshow(RED/max(max(RED)))

figure('Name', 'Red channel')
imshow(red_ch/max(max(red_ch)))

figure('Name', 'Green corrected')
imshow(GREEN/max(max(GREEN)))

figure('Name', 'Green channel')
imshow(green_ch/max(max(green_ch)))

figure('Name', 'Blue corrected')
imshow(BLUE/max(max(BLUE)))

figure('Name', 'Blue channel')
imshow(blue_ch/max(max(blue_ch)))


%Save images
red = imresize(RED,0.25);
green = imresize(GREEN,0.25);
blue = imresize(BLUE,0.25);
i = imresize(I,0.25);

filename_red = strcat(imgName(1:end-4), "red.png");
filename_green = strcat(imgName(1:end-4), "green.png");
filename_blue = strcat(imgName(1:end-4), "blue.png");
filename_i = strcat(imgName(1:end-4), "lowres.png");

imwrite(red,filename_red);
imwrite(green,filename_green);
imwrite(blue,filename_blue);



%% Function: getIntensities
function [Ir,Ig,Ib] = getIntensities(Img)
    % Extracts intensities of each color channel from an RGB image.
    % Inputs:
    %   Img - RGB image [HxWx3]
    % Outputs:
    %   Ir - intensity in the red channel
    %   Ig - intensity in the green channel
    %   Ib - intensity in the blue channel

    red_ch = double(Img(:,:,1));
    green_ch = double(Img(:,:,2));
    blue_ch = double(Img(:,:,3));

    %maximum gives much better results than average intensities -> reason:
    %mean considers pixels outside the glare points
    Ir = max(max(red_ch));
    Ig = max(max(green_ch));
    Ib = max(max(blue_ch));
    %Ir = mean(mean(red_ch));
    %Ig = mean(mean(green_ch));
    %Ib = mean(mean(blue_ch));
end
