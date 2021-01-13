% An example on how to validate the found transcription site by the eye.
% It can be runned from the workspace obtained after completed main.m

%clear all
% ----------- user inputs -----------------
%filename{1} = '.mat'; %filename, where ws of main is saved
% -----------------------------------------

%load(filename{1});
norm_img_soma = img_soma./max(img_soma, [], 'all');
norm_img_fish = img_fish./max(img_fish, [], 'all');

%% soma left
max_int_pos_left = trans_site.max_pos_left_verify;
final = [];
bb = table2array(regionprops3(bin_mask_left, 'boundingbox'));

for i = 1:size(img_soma,3)
%test1 = im2uint16(img_soma(:,:,30)./max(img_soma, [], 'all'), '16bit');
   
    test1 = norm_img_soma(:,:,i).*bin_mask_left(:,:,i)*3;
    test2 = norm_img_fish(:,:,i).*bin_mask_left(:,:,i)*6;
    final(:,:,1) = test1 ; 
    final(:,:,2) = test2;
    final(:,:,3) = zeros(size(test2));
    if abs(i-max_int_pos_left(3))<3
        
        radius=5;
        theta = linspace(0, 2*pi, round(4 * pi * radius)); % Define angles
        % Get x and y vectors for each point along the circumference.
        x = radius * cos(theta) + max_int_pos_left(1);
        y = radius * sin(theta) + max_int_pos_left(2);
            for k = 1 : length(x)
                row = round(y(k));
                col = round(x(k));
                final(row,col,3) = 255;
            end
        pause(0.5);
    end
    
    imshow(final, [], 'XData',[bb(1) bb(1)+bb(4)], 'YData',[bb(2) bb(2)+bb(5)], 'Border','tight','InitialMagnification', 500);
end

%% soma right
max_int_pos_right = trans_site.max_int_pos_right_verify;

final = [];
norm_img_soma = img_soma./max(img_soma, [], 'all');
norm_img_fish = img_fish./max(img_fish, [], 'all');
bb = table2array(regionprops3(bin_mask_right, 'boundingbox'));
for i = 1: size(img_soma,3)
%test1 = im2uint16(img_soma(:,:,30)./max(img_soma, [], 'all'), '16bit');
   
    test1 = norm_img_soma(:,:,i).*bin_mask_right(:,:,i)*3;
    test2 = norm_img_fish(:,:,i).*bin_mask_right(:,:,i)*6;
    final(:,:,1) = test1 ; 
    final(:,:,2) = test2;
    final(:,:,3) = zeros(size(test2));
    if abs(i-max_int_pos_right(3))<5
        
        radius=5;
        theta = linspace(0, 2*pi, round(4 * pi * radius)); % Define angles
        % Get x and y vectors for each point along the circumference.
        x = radius * cos(theta) + max_int_pos_right(1);
        y = radius * sin(theta) + max_int_pos_right(2);
            for k = 1 : length(x)
                row = round(y(k));
                col = round(x(k));
                final(row,col,3) = 255;
            end
        pause(0.5);
    end
    
    imshow(final, [], 'XData',[bb(1) bb(1)+bb(4)], 'YData',[bb(2) bb(2)+bb(5)], 'Border','tight','InitialMagnification', 500);

end
close all