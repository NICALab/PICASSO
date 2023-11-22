%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Three color unmixing
%%
%% Title                : PICASSO allows ultra-multiplexed fluorescence imaging of spatially overlapping proteins without reference spectra measurements 
%% Authors              : Junyoung Seo, Yeonbo Sim, Jeewon Kim, Hyunwoo Kim, In Cho, Hoyeon Nam, Young-Gyu Yoon and Jae-Byum Chang
%% Authors' Affiliation : Korea Advanced Institute of Science and Technology
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

warning('off');
addpath('./Util/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgPath = '../Data/';
filename = '3color_data.tif';

saveUnmixPath = '../Results/';
save_unmix_name = [filename(1:end-4) '_unmixed.tif'];
if (exist(saveUnmixPath, 'dir') == 0), mkdir(saveUnmixPath); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Hyper-parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
per_bg = 0; %% percentile value (0~100) of background
step_size = 0.2; %% step_size of updating alpha
maxIter = 200; %% the number of iteration
n_color = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% load input images  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_channel = single(imread([imgPath filename],'tif',1));
imgResolution = size(current_channel)
sourceIMG = zeros( imgResolution(1), imgResolution(2), n_color, 'single');

for ci=1:n_color
    current_channel = single(imread([imgPath filename],'tif',ci));
    sourceIMG(:,:,ci) = current_channel/max(current_channel(:));
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% remove background  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourceIMG_bgremoved = zeros( imgResolution(1), imgResolution(2), n_color, 'single');

for ci=1:size(sourceIMG,3)
    current_channel = sourceIMG(:,:,ci);
    estimated_background = prctile( current_channel(:), per_bg);
    current_channel = max( current_channel-estimated_background,0);
    sourceIMG_bgremoved(:,:,ci)  = current_channel;        
end    
sourceIMG_bgremoved_show = perChNorm(sourceIMG_bgremoved);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3 = figure; set(h3, 'Name', 'bgRemoved_source image_1'); imshow(sourceIMG_bgremoved_show(:,:,1), 'border','tight');
h3 = figure; set(h3, 'Name', 'bgRemoved_source image_2'); imshow(sourceIMG_bgremoved_show(:,:,2), 'border','tight');
h3 = figure; set(h3, 'Name', 'bgRemoved_source image_3'); imshow(sourceIMG_bgremoved_show(:,:,3), 'border','tight');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 3 color Unmixing %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qN = 100;
[imgDemixed, unmixing_log] = PICASSO_3C(sourceIMG_bgremoved, qN, maxIter, step_size, 0);
imgDemixed_show = perChNorm(imgDemixed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3 = figure; set(h3, 'Name', 'Unmixed image_1'); imshow( imgDemixed_show(:,:,1), 'border','tight');
h3 = figure; set(h3, 'Name', 'Unmixed image_2'); imshow( imgDemixed_show(:,:,2), 'border','tight');
h3 = figure; set(h3, 'Name', 'Unmixed image_3'); imshow( imgDemixed_show(:,:,3), 'border','tight');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% save results %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgDemixed_export = imgDemixed;
for ci=1:n_color
   current_channel = imgDemixed_export(:,:,ci);
   current_channel = current_channel/max(current_channel(:));
   imgDemixed_export(:,:,ci) = current_channel;
end

imwrite(uint16(65535*squeeze(imgDemixed_export(:,:,1))),[saveUnmixPath save_unmix_name])
for ci = 2:n_color
    imwrite(uint16(65535*squeeze(imgDemixed_export(:,:,ci))),[saveUnmixPath save_unmix_name],'WriteMode','append')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
