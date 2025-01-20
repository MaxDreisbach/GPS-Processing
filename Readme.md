# GPS-Processing
Image processing for glare-point shadowgraphy experiments (code)

by Maximilian Dreisbach (Institute of Fluid Mechanics (ISTM) - Karlsruhe Institute of Technology (KIT))


# In-situ Color Correction: P1_Extract_ISC_correction_matrix.m
This script constructs a transfer matrix for in-situ RGB correction using images captured under red, green, and blue illuminations. It calculates a correction matrix and applies it to a test image for spectral reconstruction. Corrected and original images are visualized and saved.

# Preprocess images for volumetric reconstruction: P2_Preprocess_RGB_for_vol_rec_M3_tilt_dil.m
This script processes glare-point shadowgraphy images of droplets for volumetric reconstruction. It applies spectral correction, crops and rotates images, and generates binary masks. The output includes preprocessed RGB images and binary mask files.

# Script for evaluation of droplet impingement experiments: P3_Determine_contact_angles_M7.m
This script processes high-speed images of a droplet impacting a surface. It calculates contact angles, impact velocity, and other relevant parameters to characterize the droplet dynamics. The key steps include:
- 1. Ground detection and frame of impact determination.
- 2. Ellipse fitting to the droplet contour for size and position estimation.
- 3. Computation of dynamic properties like Weber, Reynolds, Ohnesorge, and Bond numbers.
- 4. Contact angle and contact line position evaluation over time and plotting results.

# Spectral Correction and Video Creation for Shadowgraphy Images: P4_Channel_Split_and_Spectral_Correction_batch_M4_videos.m
This script performs spectral correction on glare-point shadowgraphy images and generates videos for different views (bottom, GP, and side).


If you have any questions regarding this code, please feel free to contact Maximilian (maximilian.dreisbach@kit.edu).
