%% Fractal Dimension Analysis
% Estimates box-counting fractal dimension for TIFF images.
%
% The original script limited processing to 20 images. This version
% removes the hard-coded Windows path and uses the repository data
% directory.

clc;
clear;
close all;

%% Configuration
imgDir = 'data';
imgExt = {'.tif', '.tiff'};
maxImages = 20;

%% Discover images
files = dir(imgDir);
imgFiles = {};

for i = 1:numel(files)

    if files(i).isdir
        continue;
    end

    [~, ~, ext] = fileparts(files(i).name);

    if any(strcmpi(ext, imgExt))
        imgFiles{end+1} = fullfile(imgDir, files(i).name); %#ok<SAGROW>
    end
end

if isempty(imgFiles)
    error('No TIFF images were found in %s.', imgDir);
end

if numel(imgFiles) > maxImages
    imgFiles = imgFiles(1:maxImages);
end

%% Calculate fractal dimensions
FDs = zeros(numel(imgFiles), 1);
imageNames = strings(numel(imgFiles), 1);

for k = 1:numel(imgFiles)

    I = imread(imgFiles{k});

    if size(I, 3) == 3
        I = rgb2gray(I);
    end

    BW = imbinarize(I);

    [sizes, counts] = boxSizes(BW);

    valid = counts > 0 & sizes > 0;

    if nnz(valid) < 2
        warning('Insufficient valid box-counting points for %s.', imgFiles{k});
        FDs(k) = NaN;
        continue;
    end

    coeffs = polyfit(log(1 ./ sizes(valid)), ...
                     log(counts(valid)), 1);

    FDs(k) = coeffs(1);
    [~, name, ext] = fileparts(imgFiles{k});
    imageNames(k) = string([name ext]);

    fprintf('%d (%s) FD = %.4f\n', ...
        k, imgFiles{k}, FDs(k));
end

%% Save results
resultsTable = table(imageNames, FDs, ...
    'VariableNames', {'Image','FractalDimension'});

if ~exist('results', 'dir')
    mkdir('results');
end

writetable(resultsTable, fullfile('results', 'fractal_dimensions.csv'));

disp(resultsTable);

%% Local function
function [sizes, counts] = boxSizes(BW)

    n = min(size(BW));
    maxPower = floor(log2(n));

    if maxPower < 1
        error('Image is too small for box-counting analysis.');
    end

    sizes = 2.^(maxPower:-1:1);
    counts = zeros(size(sizes));

    for i = 1:numel(sizes)

        s = sizes(i);

        nbx = ceil(size(BW,2) / s);
        nby = ceil(size(BW,1) / s);

        count = 0;

        for x = 1:nbx
            for y = 1:nby

                x0 = (x-1)*s + 1;
                x1 = min(x*s, size(BW,2));

                y0 = (y-1)*s + 1;
                y1 = min(y*s, size(BW,1));

                block = BW(y0:y1, x0:x1);

                if any(block(:))
                    count = count + 1;
                end
            end
        end

        counts(i) = count;
    end
end
