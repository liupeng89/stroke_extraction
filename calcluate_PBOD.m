clear; clc; close all;

img_path = 'com_1.jpg';

% read characters into a matrix 
IMG = imread(img_path);

% convert image to bitmap
IMG = rgb2gray(IMG);
BIT_THRESH = 127;
IMG = IMG > BIT_THRESH;

IMG = double(IMG);
nRows = size(IMG, 1); nColumns = size(IMG, 2);
index_black = find(IMG < 0.5); index_white = find(IMG > 0.5);
% 1 is black pixels and 0 is the background pixels (white).
IMG(index_black) = 1; IMG(index_white) = 0;

% figure; imagesc(IMG); colorbar;

% Control parameters
EDGE_FLAG = 2; BODY_FLAG = 1; SING_FLAG = 3; NONE_FLAG = 0;
NEARSING_FLAG = 4; CORNER_FLAG = 5; BADCORNER_FLAG = 6;

% calculate the PBOD of each pixel
index = find(IMG > 0.5);
% find black pixels coordinates in IMG
[index_rows, index_cols] = ind2sub(size(IMG), index);
% disp(index1);
num = length(index_rows);
% disp(num);
nbin = 30;
% degree of each black pixel
degree = zeros(num, 1);

pix2dis = ones(nbin, 1);
for ii = 1:nbin
    if ii ~= nbin / 2
        pix2dis(ii) = 1 / abs(cos(ii / nbin * pi));
    end
end

% disp(pix2dis);
% figure; plot(pix2dis);

tic % start timer

for k = 1:num
    
    id_row = index_rows(k);
    id_col = index_cols(k);

    count = zeros(nbin, 1);
    
    for i = 1:nbin
        theta = i * 2 * pi / nbin;
        
        new_row = id_row; new_col = id_col;
        new_row_sc = new_row; new_col_sc = new_col;
        
        dist = 0;
        
        while IMG(new_row, new_col) >= 1
            count(i) = count(i) + 1;
            dist = dist + 1;
            
            if theta < pi/4 || theta > pi*7/4
                new_row = round(-tan(theta) * dist + id_row);
                new_col = new_col + 1;
            elseif theta < pi*5/4 && theta > pi*3/4
                new_row = round(tan(theta) * dist + id_row);
                new_col = new_col - 1;
            elseif theta >= pi/4 && theta <= pi*3/4
                new_row = new_row - 1;
                new_col = round(cot(theta) * dist + id_col);
            else
                new_row = new_row + 1;
                new_col = round(-cot(theta) * dist + id_col);
            end
            
            % bound control
            if new_row > nRows || new_row < 1 || new_col > nColumns || new_col < 1
                break;
            end
        end
    end
    
    % remove small noise of peaks
    count(count < max(1/3*max(count), mean(count))) = 0;
    
    % find peaks count 
    [peaks loc] = findpeaks(count);
    numfalsepeak = 0;
    
    if length(loc) >= 2
        dist_ = zeros(length(loc)-1, 1);
        for ii = 1:length(loc)-1
            dist_(ii) = abs(loc(ii) - loc(ii+1));
        end
        dist_min = round(45 / (360 / nbin));
        falsepeak = find(dist_ < dist_min);
        numfalsepeak = length(falsepeak);
    end
    % number of peaks
    numpeaks = length(findpeaks(count)) - numfalsepeak;
    
    [peaks2 loc2] = findpeaks([count; count]);
    numfalsepeak = 0;
    if length(loc2) >= 2
        dist_ = zeros(length(loc2)-1, 1);
        for ii = 1: length(loc2)-1
            dist_(ii) = abs(loc2(ii) - loc2(ii+1));
        end
        dist_min = round(45 / (360 / nbin));
        falsepeak = find(dist_ < dist_min);
        numfalsepeak = length(falsepeak);
    end
    
    numpeaks2 = length(findpeaks([count; count])) - numfalsepeak;
    
    if numpeaks2 > 2*numpeaks
        degree(k) = numpeaks + 1;
    else
        degree(k) = numpeaks;
    end
end
toc

edge_indexs = sub2ind(size(IMG), index_rows(degree == 1), index_cols(degree == 1));
body_indexs = sub2ind(size(IMG), index_rows(degree == 2), index_cols(degree == 2));
sing_indexs = sub2ind(size(IMG), index_rows(degree >= 3), index_cols(degree >= 3));

IMG(edge_indexs) = 2;
IMG(sing_indexs) = 3;

figure; imagesc(IMG); colorbar;










