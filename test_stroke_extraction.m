clc; clear;
img_path = 'D:\DATA\kai_images\0_4326.jpg';

I = imread(img_path);
% figure; imshow(I);

% grayscale to bitmap
bit_thresh = 127;
BW = I > bit_thresh;
figure; imshow(BW);
% inverse the coler

% Get the connect components of characters.
connect_components = Components_split(BW, 8);
disp(length(connect_components));

% for i=1: length(connect_components)
%     figure; imshow(connect_components{i});
% end

% Processing each single components of character.

com_1 = connect_components{1};
figure; imshow(com_1);

com_2 = connect_components{2};
figure; imshow(com_2);
% imwrite(com_2, 'com2.png');

% get edges
edge_1 = edge(com_1, 'Canny');
figure; imshow(edge_1);

