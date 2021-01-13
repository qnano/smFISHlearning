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

% Adjusted by Pieter van Velde, June 2020.

clear all

% Add to helper functions to path
addpath(genpath('helperfunctions'))

% ------------------ User inputs ---------------------
fileName = './data/example_set_calyx.tif';     % Location dataset which needs to be analyzed.
PSFSigma=1.39;              % Nyquist size of PSF sigma in pixels.
pH1Min = 0.95;              % 1 - False detection probability.
comp_reduction = 1;         % Computational reduction in detection by applying a mask (see LLRMap3D.m)
% -----------------------------------------------------

% Extract FISH data in proper format from data.
a = bfopen(fileName);
img_FISH = double(cell2mat(permute(a{1}(:,1),[3 2 1])));

% Create mask, by selecting contour of the 'dark region' 3 times. 
img_calyx = img_FISH;
disp('Select calyx first time, right-click create mask after making polygon');
dipshow(img_calyx(:,:,round(size(img_calyx,3)/2)-7))
mask1 = roipoly;
close all
disp('Select calyx second time, right-click create mask after making polygon');
dipshow(img_calyx(:,:,round(size(img_calyx,3)/2)))
mask2 = roipoly;
close all
disp('Select calyx third time, right-click create mask after making polygon');
dipshow(img_calyx(:,:,round(size(img_calyx,3)/2)+7))
mask3 = roipoly;
close all
mask = zeros(size(img_calyx));

% Create mask calyx.
mask(:, :, [round(size(img_calyx,3)/2)-7-3:round(size(img_calyx,3)/2)-7+3])= repmat(mask1, [1,1,7]);
mask(:,:, [round(size(img_calyx,3)/2)-3:round(size(img_calyx,3)/2)+3]) = repmat(mask2, [1,1,7]);
mask(:, :, [round(size(img_calyx,3)/2)+7-3:round(size(img_calyx,3)/2)+7+3]) =repmat(mask3, [1,1,7]);


%% GLRT Based Detection, 
%with mask in there already, so only detection within the calyx are found

[coordsCam1,detParCam1,cutProcessCam1] = LLRMap3D(img_FISH,[PSFSigma 2.5*PSFSigma],[],comp_reduction,0.05,20,true,1000000, mask);

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
results.threeDPos = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1:2) coodsUnCut(maskPreFiltCam1&maskFilt1,3)];
% Lists of brightness, background and size detections
results.signal = [rawFitResultsCam1.Photons(maskPreFiltCam1&maskFilt1)];
results.bg = [rawFitResultsCam1.Bg(maskPreFiltCam1&maskFilt1)];
results.width = [rawFitResultsCam1.Sigma(maskPreFiltCam1&maskFilt1, 1)];
