clear all; close all; clc;

%% Initialize
angleSet = 0:5:175;
spatialSet = 0:1:40; % offset 
d1_conversion = 1;
d2_conversion = 1.^2;

% Input sample -- common file name
fname = 'CR-2';
I = imread([fname '_hTnI_488.tif']); % read hTnI channel
I = rgb2gray(I);
I = im2double(I);
M = imread([fname '_ManualSegment01_partial.tif']); % locate mask image
M = M > 0;

warning('off','all');
m = morph_texture_function(I,M,angleSet,spatialSet,d1_conversion,d2_conversion,fname);

plotCorrelation(m,d1_conversion,1)

% optional to create a bounding box image of the region of interest
%imwrite(m.AACT_BoundingBox,[fname '_ManualSegment_PartialCell_BB.tif']);