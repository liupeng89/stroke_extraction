% Character stroke extraction algorithm
clear; clc; close all;

% Load image
img_path = 'D:\DATA\kai_images\0_0071.jpg';

% Read image into matrix
IMG = imread(img_path);
bit_thresh = 127;
% [0, 255] to [0, 1]
IMG = IMG > bit_thresh;

% figure; imshow(IMG);

IMG = double(IMG);
n_rows = size(IMG, 1);
n_columns = size(IMG, 2);

black_index = find(IMG < 0.5);
white_index = find(IMG > 0.5);

% Change black pixel value to 1
IMG(black_index) = 1; IMG(white_index) = 0;
% figure; imshow(IMG);

% PBOD of each pixel
indexs = find(IMG > 0.5);
[index_row, index_col] = ind2sub(size(IMG), indexs);
% disp([index_row index_col])
num = length(index_row);
nbin = 30;

degree = zeros(num, 1);

tic

for k = 1:num
    row = index_row(k);
    col = index_col(k);
    
    count = zeros(nbin, 1);
    for i = 1:nbin
        theta = i*2*pi/nbin;
        
        
    end
end


toc

