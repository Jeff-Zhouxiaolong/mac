clc;clear all;
% Sample data
% texture = load('texture.mat');
% x = texture.f;
x=rgb2gray(imread(('9.png')));
[ imf_matrix ] = bemd( x );
analytic_img = hilbert2( imf_matrix(:,:,4));
analytic_img=fliplr(flipud(analytic_img));
figure
subplot 321
imagesc(x)
colormap gray;
axis equal;
axis off;
subplot 322
imagesc(imf_matrix(:,:,1))
colormap gray;
axis equal;
axis off;
subplot 323
imagesc(imf_matrix(:,:,2))
colormap gray;
axis equal;
axis off;
subplot 324
imagesc(imf_matrix(:,:,3))
colormap gray;
axis equal;
axis off;
subplot 325
imagesc(imf_matrix(:,:,4))
colormap gray;
axis equal;
axis off;
subplot 326
imagesc(analytic_img)
% imagesc((imf_matrix(:,:,1)+imf_matrix(:,:,2)+imf_matrix(:,:,3)))
colormap gray;
axis equal;
axis off;