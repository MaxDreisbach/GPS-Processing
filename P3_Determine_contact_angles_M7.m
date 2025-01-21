clear all
clc
close all

%% Script for evaluation of droplet impingement experiments
% This script processes high-speed images of a droplet impacting a surface.
% It calculates contact angles, impact velocity, and other relevant parameters
% to characterize the droplet dynamics. The key steps include:
% 1. Ground detection and frame of impact determination.
% 2. Ellipse fitting to the droplet contour for size and position estimation.
% 3. Computation of dynamic properties like Weber, Reynolds, Ohnesorge, and Bond numbers.
% 4. Contact angle and contact line position evaluation over time.
% 5. Plotting of the results
% 
% Dependencies:
% - Requires custom function: EllipseFitByTaubin
% Usage:
% - Set the `source_dir` variable to the folder containing the image sequence.
% - Adjust parameters such as `FPS`, `calibration_px_to_mm`, and `INT_CORR` as needed.
% - Enable flags `Plotting` and `Plotting2` for debugging visualizations.


%% Configuration Parameters
%Experimental conditions
hor_oscillation = true; % Set to true if the droplet exhibits a larger 
% horizontal diameter (D_H) than vertical (D_V) -> TODO: automate
calibration_px_to_mm = 87.8; % Digital image magnification [px/mm]
Rec_factor = 256; % Calibration factor for volumetric reconstruction 
% [3D-domain: 2x2x2, image domain 512px x 512px]
FPS = 7500; % frame rate of the recording

% Fluid properties for water (dynamic viscosity, density, and surface tension)
t_c = 25; %temperature [°C]
[mu_D,rho_D,sigma_D] = getFluidProperties(t_c);

%Image Processing
T_STEPS = 10; % Time steps before impact for velocity calculation; 
% adjust for balance between accuracy and frame acceleration,
% more: higher accuracy, less: lower error due to acceleration between frames,
% check bedug plots for linearity
POLYNOM_ORDER = 2; % Polynomial order for subpixel edge fitting
KERNEL_SIZE = 2; % Kernel size for subpixel fitting
OFFSET = 12; % Offset to avoid shadows or ground structure artifacts in contour detection
INT_CORR = 4; % Intensity correction factor for the blue channel
G_thresh = 0.1; % Threshold value for edge detection
EVAL_ROWS = 10; % Number of rows for evaluating contact angles

% Flags for optional plotting during debugging
Plotting = false; % Enable or disable general plots
Plotting2 = false; % Enable or disable detailed ellipse fitting plots

%% Source: image path, image name
source_dir = '../120123_SMG4_Locaton 3_Droplet 1_C001H001S0001/';
out_dir = 'eval_csv/';
out_dir_plots = 'plots/';
images = dir([source_dir,'*.tif']);  % Vector of all image filenames

%% Read first image and split cannels of RGB-image
imgName = images(1).name;   
I = imread([source_dir,imgName]);
I = double(I)/2^12;
B1 = I(:,:,3)*INT_CORR;

%% Ground detection
[BW,threshOut,Gx,Gy] = edge(B1,"Sobel",G_thresh);
V_mean = mean(Gy,2);
[m,Ground] = max(V_mean)

%Plotting
if Plotting
    imshow(B1)
    hold on
    p1 = [Ground,1];
    p2 = [Ground,1280];
    plot([p1(2),p2(2)],[p1(1),p2(1)],'Color','r','LineWidth',2)
end


%% Determine frame of impact
for i = 1:length(images)
    imgName = images(i).name;   
    I = imread([source_dir,imgName]);
    I = double(I)/2^12;
    %5px off ground to compensate tilted ground
    B = I(Ground-10:Ground,:,3)*INT_CORR;
    %B = I(:,:,3)*INT_CORR;
    %background subtraction
    B = B1(Ground-10:Ground,:) - B;
    %B = B1 - B;

    Int(i) = mean(mean(B));

    % Impact frame identified when intensity exceeds threshold
    % sometimes requires tuning
    if i>1 && Int(i) > 0.05
        f_impact = i
        break
    end 
end

if Plotting
    figure
    imshow(B)
    figure
    plot(Int)
end



%% Part 1: Determine impact velocity, droplet diameter, and droplet volume
for i = f_impact-T_STEPS:f_impact
    imgName = images(i).name;   
    I = imread([source_dir,imgName]);
    I = double(I)/2^12;
    B = I(1:Ground-OFFSET,:,3)*INT_CORR; % Region of interest above ground
    B_img = B;
    %B = B1(1:Ground-OFFSET,:) - B;
    
    %optional: median filter to remove noise
    %B = medfilt2(B);

    %image closing to remove p=1 GP
    se = strel('disk',10);
    closeBW = imclose(imcomplement(B),se);
    
    [BW,threshOut,Gx,Gy] = edge(closeBW,"Sobel",G_thresh);
    
    %% Fit ellipse to edge points
    [row_idx, col_idx] = find(BW);
    XY = [col_idx row_idx];
    
    Ax_1 = EllipseFitByTaubin(XY); 
    % fit proposed by G. Taubin (1991) - doi.org/10.1109/34.103273
    a = Ax_1(1); b = Ax_1(2)/2; c = Ax_1(3); d = Ax_1(4)/2; f = Ax_1(5)/2; g = Ax_1(6);

    center(1) = (c*d - b*f)/(b^2-a*c);
    center(2) = (a*f - b*d)/(b^2-a*c);
    sem_ax_1 = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
    sem_ax_2 = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c))));
    
    if sem_ax_1 > sem_ax_2
        sem_ax_big = sem_ax_1;
        sem_ax_small = sem_ax_2;
    else
        sem_ax_big = sem_ax_2;
        sem_ax_small = sem_ax_1;
    end

    % Get orientation angle of fitted ellipse
    if b == 0 && a < c
    phi = 0;
    elseif b == 0 && a > c
    phi = 0.5*pi;
    elseif b ~= 0 && a < c
    phi = 0.5* acot((a-c)/(2*b));
    else
    phi = 0.5*pi + 0.5* acot((a-c)/(2*b));
    end
    
    % Get measured quantities
    H(i) = center(2);
    Xpos(i) = center(1);
    Ax_1 = sem_ax_big * 2;
    Ax_2 = sem_ax_small * 2;
    Aspect(i) = Ax_1/Ax_2; % Aspect ratio of ellipse
    
    %The equivalent diameter should always be calculated as D_h^2*D_v
    %(see Zhang et al. (2019) - doi.org/10.1007/s00348-019-2712-7)
    if hor_oscillation
        D_eq(i) = (Ax_1^2 * Ax_2)^(1/3); %equivalent spherical diameter  
        Vol(i) = (Ax_1/2)^2 * (Ax_2/2) * 4/3 * pi; %volume of elipsoid (spheroid=Rotationselipsoid) with r_1 = r_2 and r_3 < r_1
    else
        D_eq(i) = (Ax_1 * Ax_2^2)^(1/3); %equivalent spherical diameter  
        Vol(i) = (Ax_1/2) * (Ax_2/2)^2 * 4/3 * pi;
    end
    
    AX_1(i) = Ax_1;
    AX_2(i) = Ax_2;
    
    % Plotting for debug
    if Plotting2
        BS = zeros(size(I,1),size(I,2));
        BS(1:size(BW,1),:) = abs(Gx);
        BW2 = double(BS);
        BW3 = cat(3,BW2,BW2,BW2);
        IC = I*10 + BW3;
        figure 
        imshow(B_img)
        hold on
        yline(center(2),'LineWidth',1,'Color',[1,0,0])
        xline(center(1),'LineWidth',1,'Color',[1,0,0])
        plot(center(1), center(2),'-o','LineWidth',3)
        viscircles(center,D_eq(i)/2,'LineWidth',1,'Color',[0,1,1])
        viscircles(center,AX_1(i)/2,'LineWidth',1,'Color',[1,0,0])
        viscircles(center,AX_2(i)/2,'LineWidth',1,'Color',[0,0,1])
    end
end

%% Save values for all considered frames
H = H(end-T_STEPS:end-1);
centerX = mean(Xpos(end-T_STEPS:end-1));
AX_1 = AX_1(end-T_STEPS:end-1);
AX_2 = AX_2(end-T_STEPS:end-1);
Vol = Vol(end-T_STEPS:end-1);
D_eq = D_eq(end-T_STEPS:end-1);

% Averaging
D_px = mean(D_eq);
D_px_std = std(D_eq);
D_phys = D_px / calibration_px_to_mm / 1000; % [m]
D_3D = D_px / Rec_factor;
Vol_px = mean(Vol);
Vol_px_std = std(Vol);
Vol_phys = Vol_px / calibration_px_to_mm^3; % [mm^3 = 10^-6 l = 1 microlitre]
Vol_3D = Vol_px / Rec_factor^3;


%Linear fit to droplets center of mass for velocity determination
x_fit = linspace(1,length(H),length(H));
P = polyfit(x_fit,H,1);
Hfit = P(1)*x_fit+P(2);
Delta_H = P(1); %difference between two frames
Delta_t = 1/FPS;
U_0 = (Delta_H / Delta_t) / calibration_px_to_mm / 1000; % [m/s]

if Plotting
    figure
    plot(H)
    hold on;
    plot(x_fit,Hfit,'r-.');
    title("Falling drop center position and linear fit")
end

%% Calculate dimensionless numbers
[We,Re,Oh,Bo] = getDimensionlessNumbers(D_phys,U_0,mu_D,rho_D,sigma_D);

% write meta-file
csvname = [out_dir,imgName(1:end-24),'_meta.csv'];
mkdir(out_dir)
varNames = ["We","Re","Oh","Bo","D","U","mu","rho","sigma"];
Table_meta = table(We,Re,Oh,Bo,D_phys,U_0,mu_D,rho_D,sigma_D,'VariableNames',varNames);
writetable(Table_meta,csvname,'Delimiter','comma')




%% Part 2: Droplet boundary detection -> contact angles & contact line position
% This section calculates the contact angles and contact line positions for a droplet
% after impact. The results include left, right, and average contact angles as well
% as the spreading factor and dynamic contact angle hysteresis.

x_lin = linspace(1,5,5);
x_lin2 = linspace(1,5,500);

for i = f_impact:length(images)
    ind = i-f_impact+1;
    T(ind) = ind/FPS;
    imgName = images(i).name;   
    I = imread([source_dir,imgName]);
    I = double(I)/2^12;
    B = I(:,:,3)*INT_CORR;
    B = imcomplement(imfill(imcomplement(B))); % fill holes in inverted image to get rod of p=1 glare point in the middle of the image
    %figure
    %imshow(B)
    l_bound = Ground-OFFSET+1;
    u_bound = l_bound-eval_rows;
    B = B(u_bound:l_bound,1:end);
    [BW,threshOut,Gx,Gy] = edge(B,"Sobel",G_thresh);
        
    %% Subpixel edge locations
    end_left = centerX - 0; % value >0 spare out glare point in the middle
    start_right = centerX + 0;
    
    for j=1:eval_rows
        left = abs(Gx(j,1:end_left));
        right = abs(Gx(j,start_right:end));
        l_ind(j) = getSubpixelEdge(left,POLYNOM_ORDER,KERNEL_SIZE);
        r_ind(j) = getSubpixelEdge(right,POLYNOM_ORDER,KERNEL_SIZE) + start_right;
    end
    
    %% linear fit to droplet boundary
    x_lin3 = linspace(1,10,10);
    %left
    c = polyfit(x_lin3,l_ind,1);
    b_left = polyval(c,x_lin3);
    CA_l(ind) = atan((b_left(end)-b_left(1))/(x_lin3(end)-x_lin3(1))) * 180/pi + 90;
    %contact line from lowest position (better: extent fit to Ground)
    CL_l(ind) = b_left(end);

    c = polyfit(x_lin3,r_ind,1);
    b_right = polyval(c,x_lin3);
    CA_r(ind) = -atan((b_right(end)-b_right(1))/(x_lin3(end)-x_lin3(1))) * 180/pi + 90;
    CL_r(ind) = b_right(end);
    
    %wetted "diameter" (to stay consistent with defition of spreading factor) - half-axis for elliptical droplet
    D(ind) = CL_r(ind) - CL_l(ind);
    
    
    if Plotting && mod(i,10)==0
        figure
        imshow(B)
        hold on
        plot([b_left(1),b_left(end)],[1,10],'Color','r','LineWidth',2)
        plot([b_right(1),b_right(end)],[1,10],'Color','r','LineWidth',2)
    end
    i
end
    
%Calculations
CA_m = 1/2*(CA_l+CA_r);
D_spread = (CL_r - CL_l) / D_px;

%Calculate contact line velocity -> Droplet capillary number -> Plot
%contact angle hysteresis (CA over capillary number)
for i = 2:length(CL_l)
    Delta_CL_l(i) = -(CL_l(i) - CL_l(i-1));
    Delta_CL_r(i) = CL_r(i) - CL_r(i-1);
    Delta_CL(i) = (Delta_CL_l(i) + Delta_CL_r(i)) / 2;
    U_CL(i) = (Delta_CL(i)/ Delta_t) / calibration_px_to_mm / 1000; % [m/s]
    Capillary_number(i) = U_CL(i) *  mu_D / sigma_D;
end


% Plot contact angle and contact line position
mkdir([out_dir_plots,imgName(1:end-24)])
% Plot contact angles
figure
plot(T,CA_l)
hold on
plot(T,CA_r)
plot(T,CA_m,'LineWidth',2)
xlabel('T [s]')
ylabel('\theta [°]')
legend('CA left','CA right','CA average')
title("Contact angles")

filename = [out_dir_plots,imgName(1:end-24),'/',imgName(1:end-24),'_contact_angles'];
saveas(gcf,filename,'pdf')

% Plot contact line position
figure
plot(T,CL_l)
hold on
plot(T,CL_r)
xlabel('T [s]')
ylabel('X [px]')
legend('CL left','CL right')
title("Contact line position")

filename = [out_dir_plots,imgName(1:end-24),'/',imgName(1:end-24),'_contact_line_position'];
saveas(gcf,filename,'pdf')

% Plot contact line position
figure
plot(T,D_spread)
xlabel('T [s]')
ylabel('D* [-]')
title("Spreading factor")

filename = [out_dir_plots,imgName(1:end-24),'/',imgName(1:end-24),'_spreading_factor'];
saveas(gcf,filename,'pdf')

% Plot contact angle hysteresis
figure 
scatter(Capillary_number,CA_m)
xlim([-0.01 0.01])
title("Dynamic contact angles")
filename = [out_dir_plots,imgName(1:end-24),'/',imgName(1:end-24),'_hoffmann_graph']
saveas(gcf,filename,'pdf')
    

csvname = [out_dir,imgName(1:end-24),'_dynamics.csv'];
mkdir(out_dir)
varNames = ["time","spreading factor",'CA left','CA right','CA average','CL left','CL right'];
Table_Spread_time = table(T',D_spread',CA_l',CA_r',CA_m',CL_l',CL_r','VariableNames',varNames);
writetable(Table_Spread_time,csvname,'Delimiter','comma')
    
    
%% Function: Subpixel Edge Location
function subPix = getSubpixelEdge(row,pol_order,kernel_size)
    % Determines subpixel-accurate edge location in an image row
    % Inputs:
    %   row - Gradient values from the image row where the edge is located
    %   pol_order - Polynomial order used for subpixel interpolation
    %   kernel_size - Size of the window around the detected peak for refinement
    % Output:
    %   subPix - Subpixel-accurate location of the edge

    x_lin = linspace(1,5,5);
    x_lin2 = linspace(1,5,5000); % higher resolution for subpixel interpolation

    len = size(row,2);
    [m,idx] = max(row);
    
    %check if maximum index is within bounds
    if idx < 3 || idx >= len-1
        idx = 3;
    end
    
    kernel = row(idx-kernel_size:idx+kernel_size);
    
    %subpixel refinement
    p = polyfit(x_lin,kernel,pol_order); % fitting polynomian
    fit_f = polyval(p,x_lin2); %subpixel interpolation
    [m,idx_subpix] = max(fit_f); %determination of subpixeel maximum
    subpix_pos  = x_lin2(idx_subpix) - 3;
    subPix = idx + subpix_pos;

end


%% Function: Fluid Properties
function [mu,rho,sigma] = getFluidProperties(T)
    % Computes material properties of water at a given temperature
    % from temperature fitting functions 
    % Input:
    %   T - Temperature in degrees Celsius (0 to 100°C)
    % Output:
    %   mu - Dynamic viscosity [Pa s]
    %   rho - Density [kg/m^3]
    %   sigma - Surface tension [N/m]
    
    % Equation for dynamic viscosity of water by Cheng (2008)
    % https://doi.org/10.1021/ie071349z
    mu = 0.001*1.790*exp((-1230-T)*T/(36100+360*T));
    
    % Equation for density of water by Volk and Kaehler (2008)
    % https://doi.org/10.1007/s00348-018-2527-y
    rho = 1000.*(1-abs((T-3.98)./615).^1.71);
    
    % DIPPR106 Equation
    % http://ddbonline.ddbst.com/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe
    T_K = T + 273.15 ; % temperature in Kelvin 
    A = 134.15; B = 1.6146; C = -2.035; D = 1.5598; T_c = 647.3; % required parameters 
    T_r = T_K/T_c;
    sigma = A *((1-T_r)^(B + C*T_r + D*T_r^2)) / 1000;
end


%% Function: Dimensionless Numbers
function [We,Re,Oh,Bo] = getDimensionlessNumbers(D,U,mu,rho,sigma)
    % Computes Weber, Reynolds, Ohnesorge, and Bond numbers
    % Inputs:
    %   D - Equivalent droplet diameter before impact [m]
    %   U - Impact velocity of the droplet [m/s]
    %   mu - Dynamic viscosity of the droplet fluid [Pa s]
    %   rho - Density of the droplet fluid [kg/m^3]
    %   sigma - Surface tension of the droplet fluid [N/m]
    % Outputs:
    %   We - Weber number, representing inertial to surface tension forces
    %   Re - Reynolds number, representing inertial to viscous forces
    %   Oh - Ohnesorge number, a dimensionless ratio of viscous to inertial and surface forces
    %   Bo - Bond number, representing gravitational to surface tension forces
    g = 9.81;
    We = D/2 * U^2 * rho / sigma;
    Re = D/2 * U * rho / mu;
    Oh = sqrt(We) / Re;
    Bo = rho * g * (D/2)^2 / sigma;
end
