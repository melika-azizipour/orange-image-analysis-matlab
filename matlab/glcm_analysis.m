%% GLCM Texture Analysis
% Computes GLCM texture features for a grayscale image and its
% SSR/MSR-enhanced versions.
%
% Features:
%   Contrast, Correlation, Energy, Homogeneity
%
% The image is divided into a 3 x 3 grid and minimum/maximum values
% and their block locations are reported.

clc;
clear;
close all;

%% Configuration
imageFile = fullfile('data', 'o1.JPG');
numBlocks = 3;
ssrSigma = 50;
msrSigmas = [15, 80, 250];

%% Read image
assert(isfile(imageFile), 'Image not found: %s', imageFile);

img = imread(imageFile);

if size(img, 3) == 3
    imgGray = rgb2gray(img);
else
    imgGray = img;
end

imgGray = im2double(imgGray);

%% Retinex preprocessing
R_ssr = ssr(imgGray, ssrSigma);
R_msr = msr(imgGray, msrSigmas);

%% Analyze blocks
analyze_blocks(imgGray, 'Gray', numBlocks, 'results');
analyze_blocks(R_ssr, 'SSR', numBlocks, 'results');
analyze_blocks(R_msr, 'MSR', numBlocks, 'results');

fprintf('GLCM analysis completed.\n');

%% Local functions
function result = ssr(img, sigma)
    I_blur = imgaussfilt(img, sigma);
    result = log(1 + img) - log(1 + I_blur);
    result = mat2gray(result);
end

function result = msr(img, sigmas)
    R = zeros(size(img));

    for i = 1:numel(sigmas)
        I_blur = imgaussfilt(img, sigmas(i));
        R = R + log(1 + img) - log(1 + I_blur);
    end

    result = R / numel(sigmas);
    result = mat2gray(result);
end

function analyze_blocks(img, label, numBlocks, outputDir)

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    [rows, cols] = size(img);

    blockH = floor(rows / numBlocks);
    blockW = floor(cols / numBlocks);

    features = {'Contrast', 'Correlation', 'Energy', 'Homogeneity'};

    maxVals = struct();
    minVals = struct();
    maxLoc = struct();
    minLoc = struct();

    for f = 1:numel(features)
        maxVals.(features{f}) = -inf;
        minVals.(features{f}) = inf;
    end

    for r = 1:numBlocks
        for c = 1:numBlocks

            rStart = (r-1) * blockH + 1;
            rEnd = r * blockH;

            cStart = (c-1) * blockW + 1;
            cEnd = c * blockW;

            block = img(rStart:rEnd, cStart:cEnd);

            glcm = graycomatrix(im2uint8(block), ...
                'NumLevels', 8, ...
                'GrayLimits', [], ...
                'Symmetric', true);

            stats = graycoprops(glcm, features);

            for f = 1:numel(features)
                val = stats.(features{f});

                if val > maxVals.(features{f})
                    maxVals.(features{f}) = val;
                    maxLoc.(features{f}) = [r, c];
                end

                if val < minVals.(features{f})
                    minVals.(features{f}) = val;
                    minLoc.(features{f}) = [r, c];
                end
            end
        end
    end

    fprintf('\n--- %s ---\n', label);

    resultRows = cell(numel(features) * 2, 5);
    row = 1;

    for f = 1:numel(features)

        fprintf('Max %s: Block (%d,%d): %.4f\n', ...
            features{f}, maxLoc.(features{f})(1), ...
            maxLoc.(features{f})(2), maxVals.(features{f}));

        fprintf('Min %s: Block (%d,%d): %.4f\n', ...
            features{f}, minLoc.(features{f})(1), ...
            minLoc.(features{f})(2), minVals.(features{f}));

        resultRows(row,:) = {label, features{f}, 'Max', ...
            maxLoc.(features{f})(1), maxVals.(features{f})};
        row = row + 1;

        resultRows(row,:) = {label, features{f}, 'Min', ...
            minLoc.(features{f})(1), minVals.(features{f})};
        row = row + 1;
    end

    T = cell2table(resultRows, ...
        'VariableNames', {'Method','Feature','Type','BlockRow','Value'});

    writetable(T, fullfile(outputDir, ...
        sprintf('glcm_%s_summary.csv', lower(label))));
end
