% Code for the publication:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015 accepted. mbc.E15-06-0448
%
% Copyright 2017, University of Oxford, United Kingdom
%
% Carlas Smith, June 2017
%
% modified by Pieter van Velde for specific purpose, April 2020

% This code requires the toolbox DIPimage and CUDA 10.2 (or higher) to be installed
% The windows watchdog that will stop the GPU execution if not disabled! 
% (TdrLevel=0,TdrDelay=180; https://msdn.microsoft.com/en-us/library/windows/hardware/ff569918(v=vs.85).aspx) 
% (http://www.diplib.org/download)

close all
clear all
gpuDevice
addpath(genpath('helperfunctions'))

% --------------------------- User inputs ----------------------------
PSFSigma= 0.73;         % nyquist size of PSF sigma in pixels
savemode = true;        % save set
fileindex = 1;          % choose file to analyze
numberofchannels = 4;   % number of imaging channels
pH1Min =0.95;           % probablity of false detection to 0.05
clusterSizeMax = inf;   % to avoid not detecting trancription side
                        % other detection parameters are default


% filenames data to be analyzed
filename{1} = '../data/example_set_transcription_focus.tif';
%filename{2} = 

% z range in which the soma exist (visually determined)(if only one soma is
% present range is [])
% format: { filename{1}: [z-range_left_soma]  [z-range_right_soma];
%           filename{2}: [z-range_left_soma]  [z-range_right_soma];
%            .
%            .
%            .
%           filename{end} [z-range_left_soma]  [z-range_right_soma]}           
z_range =   {   
                [13 49] [13 51];
            };
        z_range =   {   
                [13 49] [13 51];
            };
% --------------------------------------------------------------------

% Extract relevant channels in proper format from data.
filename=filename{fileindex};
a = bfopen(filename);
temp_img = a{1}(:,1);
alter = 1:numberofchannels:length(temp_img);
for i= 1:length(temp_img)/numberofchannels
    img_fish(:,:,i) = double(temp_img{alter(i)});
    img_soma(:,:,i) = double(temp_img{alter(i)+2});
end
data670 = img_fish;


%% GLRT Based Detection

[coordsCam1,detParCam1,cutProcessCam1] = LLRMap3D(data670,[PSFSigma],0,3);

%% Pre-Filter Detection Clusters Cam1
paramsPreFilterFits = getDefaultParamsPreFilterFits;
% paramsPreFilterFits = 
%     circularityMax: 3
%     circularityMin: 0.5000
%             pH1Max: 1
%             pH1Min: 0
%       minPixelDist: 0
%     clusterSizeMin: 0
%     clusterSizeMax: 50

% Set probablity of false detection to 0.05
paramsPreFilterFits.pH1Min = pH1Min;
paramsPreFilterFits.clusterSizeMax = clusterSizeMax; 
[ maskPreFiltCam1 ] =  preFilterFits(coordsCam1,detParCam1,paramsPreFilterFits);

% Plot pre-filtered spots
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [coordsCam1(maskPreFiltCam1,2) coordsCam1(maskPreFiltCam1,1) floor(coordsCam1(maskPreFiltCam1,3))];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 

index = double(floor(255.*(1-(1-detParCam1.pH1(maskPreFiltCam1))./(1-paramsPreFilterFits.pH1Min)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(cutProcessCam1)'; %repmat(max(cutProcessCam1,[],3),[1 1 10]))';

dipSubLoc2D(optionsLLR);
drawnow;

%% MLE Fit Intensities Cam1
paramsFit = getDefaultParamsFit;
% paramsFit       
%     FitSigma: 0
%     Iterations: 10
%     MaxCudaFits: 100000
%     PSFSigma: 1.3900
%     BoxSize: 11.3400
coodsUnCut=round(coordsCam1+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(coordsCam1,1),1) ones(size(coordsCam1,1),1) zeros(size(coordsCam1,1),1)]);
paramsFit.FitSigma = true;
paramsFit.PSFSigma = PSFSigma;
paramsFit.BoxSize = (2 * 3*PSFSigma + 1);
% Fit detections
[ rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(data670)),[coodsUnCut(:,2) coodsUnCut(:,1) coodsUnCut(:,3)],paramsFit);

% Display 3D matrix for which each 2D slice is fitted with Apixelated-gaussian
dipshow(rawFitResultsCam1.ROIStack,'lin')
drawnow;

%%
paramsFilterFits = getDefaultParamsFilterFits;
% paramsFilterFits
%       MinCRLBSTD: 0
%       MinPhotons: 0
%            MinBg: 0
%       MaxCRLBSTD: Inf
%        MaxPFAValue: 0.05
%       MaxPhotons: Inf
%            MaxBg: Inf
%     MinPixelDist: 0

% Single Pixel Likelihood Ratio estimates can also be
% used as filtering
paramsFilterFits.MaxPFAValue=1;
[ maskFilt1 ] =  filterFits(rawFitResultsCam1,paramsFilterFits);

% Plot filtered spots
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1) rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2) rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 

index = double(floor(255.*(min(rawFitResultsCam1.PFA(maskPreFiltCam1&maskFilt1)./paramsFilterFits.MaxPFAValue,1)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(data670)';
dipSubLoc2D(optionsLLR);


%% display fit results 
figure
C=jet(256); 
index = double(floor(255.*(1-(1-detParCam1.pH1(maskPreFiltCam1&maskFilt1))./(1-paramsPreFilterFits.pH1Min)))+1);
BoxColor = C(index,:);
subplot(3,2,3:4)
scatter3(rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1),rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2),rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1),60,BoxColor,'marker','.')
colormap(jet);
caxis([paramsPreFilterFits.pH1Min 1]);
c = colorbar;
ylabel(c,'Detection Probability')
xlabel('x-position [pixel]')
ylabel('y-position [pixel]')
zlabel('t-frame [pixel]')
ntitle('Average Detection Probablity')

subplot(3,2,5:6)
index = double(floor(255.*(min(rawFitResultsCam1.PFA(maskPreFiltCam1&maskFilt1)./paramsFilterFits.MaxPFAValue,1)))+1);
BoxColor = C(index,:);
scatter3(rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1),rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2),rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1),60,BoxColor,'marker','.')
colormap(jet);
caxis([0 paramsFilterFits.MaxPFAValue]);
c = colorbar;
ylabel(c,'False Alarm Probability')
xlabel('x-position [pixel]')
ylabel('y-position [pixel]')
zlabel('t-frame [pixel]')
ntitle('Point Estimate of False Alarm Probablity')

subplot(3,2,1)
histogram(rawFitResultsCam1.Photons,10)
ylabel('Count [#]')
xlabel('Signal [#Photons]')
subplot(3,2,2)
histogram(rawFitResultsCam1.Bg,10)
ylabel('Count [#]')
xlabel('Background [#Photons/pixel]')

% List of 3D positions x,y by MLE and z by center of mass of signification
% hypothesis tested pixels
threeDPos = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1:2) coodsUnCut(maskPreFiltCam1&maskFilt1,3)];

% list of intensity, backgrount and width of spots plus its STD.
photons_pp = [rawFitResultsCam1.Photons(maskPreFiltCam1&maskFilt1)];
bg_pp = [rawFitResultsCam1.Bg(maskPreFiltCam1&maskFilt1)];
width_pp = [rawFitResultsCam1.Sigma(maskPreFiltCam1&maskFilt1, :)];
photons_pp_std = [rawFitResultsCam1.Photons_STD(maskPreFiltCam1&maskFilt1)];
bg_pp_std = [rawFitResultsCam1.Bg_STD(maskPreFiltCam1&maskFilt1)];
width_pp_std = [rawFitResultsCam1.Sigma_STD(maskPreFiltCam1&maskFilt1, :)];

%% Create mask of soma

% manually select largest area covered by the somas
close 'all'
disp('Select left soma, right-click create mask after making polygon');
dipshow(img_soma)
%dipshow(im2uint16(sum(img_soma,3)./max(sum(img_soma,3), [], 'all')), '16bit')
img_soma_left = roipoly;
img_soma_left = repmat(img_soma_left,[1,1,size(img_soma,3)]);
close 'all'
disp('Select right soma, right-click create mask after making polygon');
dipshow(img_soma)
img_soma_right = roipoly;
img_soma_right = repmat(img_soma_right,[1,1,size(img_soma,3)]);
disp('done');
close 'all'

% difference of gaussians to enhance edges, while suppresing noise
img_temp = double(smooth(img_soma,1)-smooth(img_soma,5));

% z range in array
z_stack_left = zeros(size(img_soma_left));
z_stack_left(:,:,z_range{fileindex,1}(1):z_range{fileindex,1}(2)) = repmat(ones(size(img_soma_left,1),size(img_soma_left,2)),[1,1,z_range{fileindex,1}(2)-z_range{fileindex,1}(1)+1]);

z_stack_right = zeros(size(img_soma_right));
z_stack_right(:,:,z_range{fileindex,2}(1):z_range{fileindex,2}(2)) = repmat(ones(size(img_soma_right,1),size(img_soma_right,2)),[1,1,z_range{fileindex,2}(2)-z_range{fileindex,2}(1)+1]);

% threshold image
mask_soma_left = double(img_soma_left.*img_temp.*z_stack_left > dip_mean(img_temp,logical(img_soma_left.*z_stack_left) ,[1,1,1])+dip_standarddeviation(img_temp,logical(img_soma_left.*z_stack_left) ,[1,1,1]));
mask_soma_right = double(img_soma_right.*img_temp.*z_stack_right > dip_mean(img_temp,logical(img_soma_right.*z_stack_right),[1,1,1])+dip_standarddeviation(img_temp,logical(img_soma_right.*z_stack_right),[1,1,1]));

% create mask for left cell
bin_mask_left= zeros(size(img_soma));
for i=z_range{fileindex,1}(1): z_range{fileindex,1}(2)
   
    if  nnz(logical(squeeze(mask_soma_left(:,:,i))))<=3
        % nothing, keep zero
    else
        C = convhull(convhull(dip_image(logical(squeeze(mask_soma_left(:,:,i))))));
        C = double(C);
        CP = ceil(contour(C,1));
        CP = CP(:,2:end);
        for j=1:length(CP)
           bin_mask_left(CP(2,j),CP(1,j),i) = 1;
        end
        bin_mask_left(:,:,i) = imfill(bin_mask_left(:,:,i), 'holes');
    end
end

% create mask for right cell
bin_mask_right= zeros(size(img_soma));
for i=z_range{fileindex,2}(1): z_range{fileindex,2}(2)
    if  nnz(logical(squeeze(mask_soma_right(:,:,i))))<=3
        % nothing, keep zero
    else
        C = convhull(convhull(dip_image(logical(squeeze(mask_soma_right(:,:,i))))));
        C = double(C);
        CP = ceil(contour(C,1));
        CP = CP(:,2:end);
        for j=1:length(CP)
           bin_mask_right(CP(2,j),CP(1,j),i) = 1;
        end
        bin_mask_right(:,:,i) = imfill(bin_mask_right(:,:,i), 'holes');
    end
end

% compute Properties of all Detections of the left soma
[max_int_pos_left, dif_lim.width_without_max_left, dif_lim.photons_without_max_left,...
            dif_lim.bg_without_max_left] = ...
    compute(bin_mask_left,threeDPos, photons_pp, bg_pp, width_pp);

% compute Properties of all Detections of the right soma
[max_int_pos_right, dif_lim.width_without_max_right, dif_lim.photons_without_max_right,...
            dif_lim.bg_without_max_right] = ...
    compute(bin_mask_right,threeDPos, photons_pp, bg_pp, width_pp);

%%

% Estimate sigma for the transcription site.
[psf_trans_site_left] = get_psf_transcription_site(max_int_pos_left, img_fish);
[psf_trans_site_right] = get_psf_transcription_site(max_int_pos_right, img_fish);
close all


% Redo whole detection procedure (very inefficient, better to only pass roi
% of transcription site to MLE function.), with estimated psf,
% and only get information of transcription site.
drawnow;
data670 = img_fish();
drawnow;
gpuDevice
[trans_site.max_pos_left_verify, trans_site.max_width_left, ...
    trans_site.max_photons_left, trans_site.max_bg_left] ...
    = detection_transcription(data670, psf_trans_site_left, pH1Min, clusterSizeMax, bin_mask_left);
    
[trans_site.max_int_pos_right_verify, trans_site.trans_max_int_width_right, ...
    trans_site.trans_max_photons_right,trans_site.max_bg_rigth] ...   
    = detection_transcription(data670, psf_trans_site_right, pH1Min, clusterSizeMax, bin_mask_right);


close all

