%% HSV-Based Color Correction
% Converts selected RGB image regions to HSV and applies a simple
% correction to reduce green hue influence.
%
% NOTE:
% Despite the original project label "CILAB/CIELAB", the supplied
% algorithm uses HSV (Hue, Saturation, Value), not CIELAB (L*a*b*).
% This file therefore uses an explicit HSV-based name.

clc;
clear;
close all;

%% Configuration
fileList = { ...
    'o1.JPG','o2.JPG','o3.JPG','o4.JPG', ...
    'o6.JPG','o7.JPG','o8.JPG','o9.JPG','o10.JPG','o11.JPG'};

rectSize = 200;
inputDir = 'data';
outputDir = 'results';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

allRects = [];
summary = cell(numel(fileList), 7);

%% Process images
for k = 1:numel(fileList)

    imageFile = fullfile(inputDir, fileList{k});

    if ~isfile(imageFile)
        warning('Image not found: %s. Skipping.', imageFile);
        continue;
    end

    img = imread(imageFile);

    figure(1);
    imshow(img, 'InitialMagnification', 'fit');
    title(sprintf('Orange %d - Draw ROI and press Enter', k));

    cropImg = imcrop(img);

    if isempty(cropImg)
        warning('No ROI selected for image %s.', fileList{k});
        continue;
    end

    hsvImg = rgb2hsv(cropImg);

    H = hsvImg(:,:,1);
    S = hsvImg(:,:,2);
    V = hsvImg(:,:,3);

    originalMeanH = mean2(H);
    originalMeanS = mean2(S);
    originalMeanV = mean2(V);

    % Reduce saturation in the green hue range.
    greenMask = (H > 0.22 & H < 0.45);
    S(greenMask) = S(greenMask) * 0.3;

    % Shift hue slightly.
    H = H - 0.02;
    H(H < 0) = H(H < 0) + 1;

    hsvImg(:,:,1) = H;
    hsvImg(:,:,2) = S;

    newRGB = hsv2rgb(hsvImg);

    meanR = mean2(newRGB(:,:,1));
    meanG = mean2(newRGB(:,:,2));
    meanB = mean2(newRGB(:,:,3));

    rect = zeros(rectSize, rectSize, 3);
    rect(:,:,1) = meanR;
    rect(:,:,2) = meanG;
    rect(:,:,3) = meanB;

    allRects = cat(2, allRects, rect);

    figure;
    imshow(rect);
    title(sprintf('Orange %d - Green Hue Corrected - %s', k, fileList{k}));

    summary{k,:} = {fileList{k}, originalMeanH, originalMeanS, ...
        originalMeanV, meanR, meanG, meanB};
end

%% Save summary
validRows = ~cellfun(@isempty, summary(:,1));
summaryTable = cell2table(summary(validRows,:), ...
    'VariableNames', {'Image','MeanH','MeanS','MeanV', ...
                      'CorrectedMeanR','CorrectedMeanG','CorrectedMeanB'});

writetable(summaryTable, fullfile(outputDir, 'hsv_color_summary.csv'));

if ~isempty(allRects)
    figure;
    imshow(allRects);
    title('All Orange Mean Colors After Green Hue Correction');
    exportgraphics(gcf, fullfile(outputDir, 'corrected_orange_colors.png'), ...
        'Resolution', 300);
end

disp('HSV-based color analysis completed.');
