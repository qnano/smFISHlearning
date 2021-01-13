function [max_int_pos_verify, trans_max_int_width,...
    trans_max_photons, trans_max_bg]...
    = detection_transcription(data670, PSFSigma, pH1Min, clusterSizeMax, bin_mask)
% detection_transcription: Return properties of the transcription site and
% perform detection according:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015 accepted. mbc.E15-06-0448
%
% Copyright 2017, University of Oxford, United Kingdom
%
% Carlas Smith, June 2017
%
% modified by Pieter van Velde for specific purpose, June 2020

% SYNOPSIS:
%[max_int_pos_verify, trans_max_int_width,...
%    trans_max_photons, trans_max_bg]...
%    = detection_transcription(data670, PSFSigma, pH1Min, clusterSizeMax, bin_mask
% 
% PARAMETERS:
%     data670: smFISH data
% 
%     PSFSigma: Estimated width of transcription site
% 
%     pH1Min: 1 - probabilty false positive detection
% 
%     clusterSizeMax: Max size detection cluster
% 
%     bin_mask: 3D mask of soma 
% 
% 
% OUTPUTS:
%   max_int_pos: Position transcription site
%   max_int_width: Width transcription site             
%   max_photons: Intensity transcription site
%   max_int_bg: Background transcitpion site



%% GLRT Based Detection

[coordsCam1,detParCam1,cutProcessCam1] = LLRMap3D(data670,[PSFSigma],[],3);

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


%% MLE Fit Intensities Cam1
paramsFit = getDefaultParamsFit;
% paramsFit       
%     FitSigma: 0
%     Iterations: 10
%     MaxCudaFits: 100000
%     PSFSigma: 1.3900
%     BoxSize: 11.3400
coodsUnCut=round(coordsCam1+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(coordsCam1,1),1) ones(size(coordsCam1,1),1) zeros(size(coordsCam1,1),1)]);
paramsFit.FitSigma=true;
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

% List of 3D positions x,y by MLE and z by center of mass of signification
% hypothesis tested pixels
threeDPos = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1:2) coodsUnCut(maskPreFiltCam1&maskFilt1,3)];

% list of intensity, backgrount and width of spots.
photons_pp = [rawFitResultsCam1.Photons(maskPreFiltCam1&maskFilt1)];
bg_pp = [rawFitResultsCam1.Bg(maskPreFiltCam1&maskFilt1)];
width_pp = [rawFitResultsCam1.Sigma(maskPreFiltCam1&maskFilt1, :)];


% find detections in mask
X = findcoord(dip_image(bin_mask));
[idx, dist] = knnsearch(X,threeDPos,'dist','cityblock','k',1);
points_in_mask = dist <= 1;

% find parameters of detections in mask
pos_in_mask = threeDPos(points_in_mask,:);
photons_in_mask = photons_pp(points_in_mask);
bg_in_mask = bg_pp(points_in_mask);

% width defined as magnitidue width in x and y direction
width_in_mask = width_pp(points_in_mask, :);
width_mag = sqrt(width_in_mask(:,1).^2+ width_in_mask(:,2).^2);

% find values of transcripion site
[~, index] = max((photons_in_mask));
trans_max_bg = bg_in_mask(index);
max_int_pos_verify = pos_in_mask(index,:);
trans_max_int_width = width_mag(index);
trans_max_photons = max(photons_in_mask);

end

