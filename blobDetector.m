clc
close all
clear all

inp = imread("img.jpg");

img = rgb2gray(inp);

[ keyPts keyPtsLoc DoG2 ] = getSIFTkeyPoints(img);

imshow(img);
hold on;
plot(keyPtsLoc(2,:), keyPtsLoc(1,:), 'r.');
