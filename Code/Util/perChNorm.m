%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Title                : PICASSO allows ultra-multiplexed fluorescence imaging of spatially overlapping proteins without reference spectra measurements 
%% Authors              : Junyoung Seo, Yeonbo Sim, Jeewon Kim, Hyunwoo Kim, In Cho, Hoyeon Nam, Young-Gyu Yoon and Jae-Byum Chang
%% Authors' Affiliation : Korea Advanced Institute of Science and Technology
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function img_out = perChNorm(img_in)

n_ch = size(img_in,3);
img_out = 0*img_in;

for ci=1:n_ch
   current_channel = img_in(:,:,ci);   
   normFactor = prctile(current_channel(:), 99.5);
   current_channel = current_channel/(normFactor+eps);
   img_out(:,:,ci) = current_channel;
end