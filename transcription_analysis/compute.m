function [max_int_pos, width_without_max, photons_without_max,...
            bg_without_max] = ...
            compute(bin_mask,threeDPos, photons_pp, bg_pp, width_pp)            
% compute: Extract information of all spots, excluding transcription site,
% defined as the highest intensity detection.
% SYNOPSIS:
%  [max_int_pos, width_without_max, photons_without_max,...
%            bg_without_max] = ...
%            compute(bin_mask,threeDPos, photons_pp, bg_pp, width_pp)
% 
% PARAMETERS:
%     bin_mask: 3D mask of soma
% 
%     threeDPos: Positions detections
% 
%     photons_pp: Intensity detections
% 
%     width_pp: Width detections
% 
%     bg_pp: Background detections
% 
% 
% OUTPUTS:
%   max_int_pos: Position transcription site
%   width_without_max: Width diffraction limitd spots            
%   photons_without_max: Intensity diffraction limitd spots 
%   bg_without_max: Background diffraction limitd spots 

        
% find detections in mask
X = findcoord(dip_image(bin_mask));
[~, dist] = knnsearch(X,threeDPos,'dist','cityblock','k',1);
points_in_mask = dist <= 1;

% find parameters of detections in mask
pos_in_mask = threeDPos(points_in_mask,:);
photons_in_mask = photons_pp(points_in_mask);
bg_in_mask = bg_pp(points_in_mask);

% width defined as magnitidue width in x and y direction
width_in_mask = width_pp(points_in_mask, :);
width_mag = sqrt(width_in_mask(:,1).^2+ width_in_mask(:,2).^2);

% find transcription site (now defined as detection with max number of photons)
photons_div_bg = photons_in_mask./bg_in_mask;
[~, index] = max((photons_in_mask));
max_int_pos = pos_in_mask(index,:);

% Define arrays of properties single molecules, without transcription site
photons_without_max = photons_in_mask; photons_without_max(index)=[];
bg_without_max = bg_in_mask; bg_without_max(index)=[];
width_without_max = width_mag; width_without_max(index)=[];

end

