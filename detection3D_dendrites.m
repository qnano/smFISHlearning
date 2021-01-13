% Code for the publication:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015 accepted. mbc.E15-06-0448
%
% Copyright 2017, University of Oxford, United Kingdom
%
% Carlas Smith, June 2017
%
% This code requires the toolbox DIPimage and CUDA 10.2 to be installed
% The windows watchdog that will stop the GPU execution if not disabled! 
% (TdrLevel=0,TdrDelay=180; https://msdn.microsoft.com/en-us/library/windows/hardware/ff569918(v=vs.85).aspx) 
% (http://www.diplib.org/download; https://developer.nvidia.com/cuda-toolkit-65)

% Adjusted by Pieter van Velde, August 2020.

% ------------------ User inputs ---------------------
fileName = './data/example_set_dendrites.r3d';     % Location dataset which needs to be analyzed.
PSFSigma=1.39;              % Nyquist size of PSF sigma in pixels.
numberOfChannels = 4;       % Number of channels data.
fish_channel = 1;           % Channel corresponding to the FISH data.
dendrite_channel = 3;       % Channel corresponding to dendrite data.
pH1Min = 0.95;              % 1 - False detection probability.
comp_reduction = 3;         % Computational reduction in detection by applying a mask (see LLRMap3D.m)
% ----------------------------------------------------

% Add to helper functions to path
addpath(genpath('./helperfunctions'))

% Extract relevant channels in proper format from data.
a = bfopen(fileName);
img = double(cell2mat(permute(a{1}(:,1),[3 2 1])));
img = reshape(img,[size(img,1) size(img,2) size(img,3)/numberOfChannels numberOfChannels]);
img_FISH = img(:,:,:,fish_channel);
img_Dendrites = img(:,:,:,dendrite_channel);



%% GLRT Based Detection

[coordsCam1,detParCam1,cutProcessCam1] = LLRMap3D(img_FISH,[PSFSigma 2.5*PSFSigma],[],comp_reduction);

%% Pre-Filter Detection Clusters Cam1
paramsPreFilterFits = getDefaultParamsPreFilterFits;
% paramsPreFilterFits = 
%     circularityMax: 3
%     circularityMin: 0.5000
%             pH1Max: 1
%             pH1Min: 0
%       minPixelDist: 0
%     clusterSizeMin: 0
%     clusterSizeMax: 100

% Set probablity of false detection to 0.05%
paramsPreFilterFits.pH1Min =pH1Min;
[ maskPreFiltCam1 ] =  preFilterFits(coordsCam1,detParCam1,paramsPreFilterFits);

% Plot pre-filtered spots
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [coordsCam1(maskPreFiltCam1,2) coordsCam1(maskPreFiltCam1,1) floor(coordsCam1(maskPreFiltCam1,3))];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 

index = double(floor(255.*(1-(1-detParCam1.pH1(maskPreFiltCam1))./(1-paramsPreFilterFits.pH1Min)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(cutProcessCam1)';

dipSubLoc2D(optionsLLR);
drawnow;

%% MLE Fit Intensities Cam1
paramsFit = getDefaultParamsFit;
% paramsFit       
%     FitSigma: 0
%     Iterations: 10
%     MaxCudaFits: 100000
%     PSFSigma: 1.3900
%     BoxSize: 3*2*PSFSigma+1
coodsUnCut=round(coordsCam1+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(coordsCam1,1),1) ones(size(coordsCam1,1),1) zeros(size(coordsCam1,1),1)]);
paramsFit.FitSigma=true;

% Fit detections
[ rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(img_FISH)),[coodsUnCut(:,2) coodsUnCut(:,1) coodsUnCut(:,3)],paramsFit);

% Display 3D matrix for which each 2D slice is fitted with a pixelated-gaussian.
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

% Don't use single pixel likelihood ratio estimates as filter.
paramsFilterFits.MaxPFAValue=1;

% Create filter mask.
[ maskFilt1 ] =  filterFits(rawFitResultsCam1,paramsFilterFits);

% Plot filtered spots.
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1) rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2) rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 

index = double(floor(255.*(min(rawFitResultsCam1.PFA(maskPreFiltCam1&maskFilt1)./paramsFilterFits.MaxPFAValue,1)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(img_FISH)';
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
% hypothesis tested pixels. 
threeDPos = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1:2) coodsUnCut(maskPreFiltCam1&maskFilt1,3)];
% Lists of brightness, background and size detections
signal = [rawFitResultsCam1.Photons(maskPreFiltCam1&maskFilt1)];
bg = [rawFitResultsCam1.Bg(maskPreFiltCam1&maskFilt1)];
width = [rawFitResultsCam1.Sigma(maskPreFiltCam1&maskFilt1, 1)];

%% Define mask of the dendrites.
close all
% Select edges of largest area covered by dendrites, scroll to image.
disp('Select left dendrites, right-click create mask after making polygon');
dipshow(img_Dendrites)
img_Dendrites_left = roipoly;
img_Dendrites_left = repmat(img_Dendrites_left,[1,1,size(img_Dendrites,3)]);
disp('Select right dendrites, right-click create mask after making polygon');
dipshow(img_Dendrites)
img_Dendrites_right = roipoly;
img_Dendrites_right = repmat(img_Dendrites_right,[1,1,size(img_Dendrites,3)]);
close all

% Difference of Gaussians to enhance edges, while suppresing noise.
img_temp = double(smooth(img_Dendrites,1)-smooth(img_Dendrites,5));

% Create mask dendrites.
mask_Dendrites_left = img_Dendrites_left.*img_temp > dip_mean(img_temp,img_Dendrites_left,[1,1,1])+dip_standarddeviation(img_temp,img_Dendrites_left,[1,1,1]);
mask_Dendrites_right = img_Dendrites_right.*img_temp > dip_mean(img_temp,img_Dendrites_right,[1,1,1])+dip_standarddeviation(img_temp,img_Dendrites_right,[1,1,1]);
mask = logical(mask_Dendrites_left) | logical(mask_Dendrites_right);

% find all detections inside mask
X = findcoord(dip_image(mask));
[idx, dist] = knnsearch(X,threeDPos,'dist','cityblock','k',1);
points_in_mask = dist <= 1;

% filter out only the detections inside the dendrites
results.threeDPos = threeDPos(points_in_mask,:);
% Lists of brightness, background and size detections
results.signal = signal(points_in_mask);
results.bg = bg(points_in_mask);
results.width = width(points_in_mask);

