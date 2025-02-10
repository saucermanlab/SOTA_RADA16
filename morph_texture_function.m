function metrics = morph_texture_function(aact,mask,angleSet,spatialSet,d1_conversion,d2_conversion,rootname)
%% FUNCTION DESCRIPTION
% Written by Matthew Sutcliffe doi: 10.1038/s41598-018-19539-z
% Adapted by Alice Luanpaisanon
% This function computes cell morphology features (cell area, elongation,
% circularity, eccentricity) and sarcomere organization features (sarcomere
% organization score, sarcomere length, cell-sarcomere misalignment) for
% every cell in the image.
%
% INPUTS
% hTnT : hTnI channel image for cardiomyocyte staining
% mask : segmentation image with same dimensions as aact. Each pixel is labeled as belonging to a cell
%   (1 = cell #1, 2 = cell #2, etc.) or background = 0.
% angleSet : array of orientations at which to measure Haralick
%   correlation, typically 0:1:179
% spatialSet : array of pixel offsets at which to measure Haralick
%   correlation, typically 0:1:40
% d1_conversion : pixels to microns conversion for a single dimension
% d2_conversion : pixels to microns conversion for 2 dimensions (i.e. area)
% rootname : image file root name for hTnI and mask
%
% OUTPUTS
% metrics : structure with all measured features, each of which contains an
%   array, where each element is the measured feature for the indexed cell
%% Split image
metrics = regionprops(mask,aact,'all');

for iCell = 1:numel(metrics)
    metrics(iCell).RootName = rootname;
    
    metrics(iCell).Area_um2 = metrics(iCell).Area.*d2_conversion;
    metrics(iCell).MajorAxisLength_um = metrics(iCell).MajorAxisLength.*d1_conversion;
    metrics(iCell).MinorAxisLength_um = metrics(iCell).MinorAxisLength.*d1_conversion;
    
    iMask = (mask == iCell);
    
    % Crop image to only include entire cell
    iaact = aact(:,:,1);
    iaact(~any(iMask,2),:) = [];
    iaact(:,~any(iMask,1)) = [];

    iMask(~any(iMask,2),:) = [];
    iMask(:,~any(iMask,1)) = [];

    metrics(iCell).AACT_BoundingBox = iaact;
    
    % Calculate elongation and circularity
    metrics(iCell).Elongation = metrics(iCell).MajorAxisLength./metrics(iCell).MinorAxisLength;
    metrics(iCell).Circularity = (4.*pi.*metrics(iCell).Area)./(metrics(iCell).Perimeter)^2;
    
    % Set up for sub-pixel resolution interpolation
    spatialSetInterp = 0:0.1:max(spatialSet);
    
    % Create offset matrix for co-occurrence matrix calculation
    offsets = [zeros(numel(spatialSet),1) spatialSet'; -spatialSet' zeros(numel(spatialSet),1)];
    
    % Initialize m x n matrix, where m = number of angles and n =
    % number of interpolated spatial pairs
    correlationdata = zeros(numel(angleSet),numel(spatialSetInterp));
    
    % For efficiency, calculate together angles that are 90 degrees apart
    halfangleSet = angleSet(1:(numel(angleSet)/2));
    
    % Calculate correlation for each angle and spatial distance
    angleIndex = 1;
    for iAngle = halfangleSet
        rotateI = imrotate(iaact,-iAngle);
        rotateM = imrotate(iMask,-iAngle);
        
        rotateI(~rotateM) = NaN;
        
        glcm = graycomatrix(rotateI,'offset',offsets);
        stats = graycoprops(glcm);
        stats.Correlation(isnan(stats.Correlation)) = 0;
        correlationdata(angleIndex,:) = interp1(spatialSet,stats.Correlation(1:numel(spatialSet)),spatialSetInterp,'spline');
        correlationdata(angleIndex+numel(angleSet)/2,:) = interp1(spatialSet,stats.Correlation(numel(spatialSet)+1:end),spatialSetInterp,'spline');
        angleIndex = angleIndex + 1;
        
    end
    
    metrics(iCell).CorrelationData = correlationdata;
    
    amplitude = zeros(size(correlationdata,1),1);
    lengthIndex = zeros(size(correlationdata,1),1);
% Initialize arrays to store amplitudes and spatial offsets for each angle
allAmplitudes = zeros(size(correlationdata,1), 1);
allSpatialOffsets = zeros(size(correlationdata,1), 1);

%%
% Loop over each angle
    for iAngle = 1:size(correlationdata,1)
           % Find peaks in the correlation curve
        [pks,locs,~,p] = findpeaks(correlationdata(iAngle,:));
% Check if peaks are found
        if ~isempty(pks)
% Identify the maximum peak and its location
            [pMax,locMax] = max(p);
% Store the maximum peak as amplitude for the current angle
            amplitude(iAngle) = pMax;
 % Store the corresponding spatial offset for the maximum peak
            lengthIndex(iAngle) = spatialSetInterp(locs(locMax));
% Store the amplitude and spatial offset for all angles
            allAmplitudes(iAngle) = pMax;
            allSpatialOffsets(iAngle) = spatialSetInterp(locs(locMax));
        end
    end
    % Update the metrics structure to include the new fields
    metrics(iCell).AllAmplitudes = allAmplitudes;
    metrics(iCell).AllSpatialOffsets = allSpatialOffsets;
   
    [sarcOrgScore,sarcIndex] = max(amplitude);
    
    metrics(iCell).SarcomereOrganizationScore = sarcOrgScore;
    metrics(iCell).SarcomereOrientationIndex = sarcIndex;
    metrics(iCell).SarcomereOrientation = angleSet(sarcIndex);
    metrics(iCell).CellSarcomereAlignment = min(180-mod(metrics(iCell).Orientation-metrics(iCell).SarcomereOrientation,180),...
        mod(metrics(iCell).Orientation-metrics(iCell).SarcomereOrientation,180));
    
    if lengthIndex(sarcIndex) > 0
        metrics(iCell).SarcomereLength = lengthIndex(sarcIndex);
    else
        metrics(iCell).SarcomereLength = 0;
    end
    metrics(iCell).SarcomereLength_um = metrics(iCell).SarcomereLength.*d1_conversion;
    
    if sarcIndex+floor(numel(angleSet)./2) <= numel(angleSet)
        widthIndex = sarcIndex+floor(numel(angleSet)./2);
    else
        widthIndex = sarcIndex-floor(numel(angleSet)./2);
    end
    
    if lengthIndex(widthIndex) > 0
        metrics(iCell).SarcomereWidth = lengthIndex(widthIndex);
    else
        metrics(iCell).SarcomereWidth = 0;
    end
    metrics(iCell).SarcomereWidth_um = metrics(iCell).SarcomereWidth.*d1_conversion;
    
end


metrics = orderfields(metrics);
end