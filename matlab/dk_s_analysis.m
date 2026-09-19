%% K/S and d(K/S)/dLambda Analysis
% Calculates Kubelka-Munk K/S and its wavelength derivative
% for multiple samples.
%
% Input:
%   wavelength - wavelength vector in nm
%   R          - reflectance matrix in percent (rows = wavelengths,
%                columns = samples)
%
% Output:
%   dk_s_table.csv
%   dk_s_derivative.png
%
% MATLAB Image/Signal Toolboxes are not required for this script.

clc;
clear;
close all;

%% Input data
wavelength = (400:10:750)';

samples = {'01m','02m','03m','04m','05m','06m','07m','08m','09m','010m','011m'};

R = [ ...
5.85 4.98 5.11 5.90 4.69 4.37 5.09 6.08 5.37 4.62 5.07
6.06 5.16 5.35 6.20 4.69 4.62 5.55 6.30 5.65 4.80 5.28
6.20 5.29 5.60 6.45 4.63 4.72 5.83 6.34 5.68 4.79 5.30
6.38 5.48 5.90 6.72 4.62 4.85 6.10 6.38 5.70 4.78 5.32
6.55 5.62 6.14 7.00 4.62 4.97 6.42 6.47 5.78 4.81 5.35
6.70 5.81 6.41 7.40 4.65 5.10 6.76 6.53 5.85 4.83 5.42
6.97 6.09 6.78 7.99 4.68 5.33 7.33 6.63 5.91 4.87 5.47
7.17 6.34 7.09 8.45 4.69 5.46 7.80 6.67 5.96 4.88 5.50
7.27 6.41 7.21 8.73 4.75 5.51 7.99 6.73 6.02 4.91 5.55
7.85 7.07 8.05 9.99 4.83 6.00 9.34 6.85 6.16 5.02 5.67
9.89 9.28 10.37 13.31 5.13 7.40 13.18 7.14 6.49 5.27 5.93
13.98 13.81 14.74 19.08 5.83 10.24 20.34 7.90 7.37 6.01 6.63
20.16 20.88 21.71 26.85 7.56 15.12 30.07 9.61 9.22 7.54 8.15
27.17 28.37 29.35 34.23 10.95 20.85 38.87 12.69 12.37 10.30 10.89
34.32 35.25 36.52 40.40 16.50 26.41 45.83 17.24 16.84 14.39 15.00
40.33 40.65 42.11 44.81 23.60 31.01 50.40 22.54 21.82 19.34 19.82
44.86 44.35 45.77 47.19 30.88 34.27 52.70 28.08 26.90 24.77 24.90
51.32 49.41 51.02 50.89 39.32 38.59 56.41 34.66 33.03 31.44 30.99
56.57 52.99 54.67 52.93 47.48 41.55 58.29 41.07 39.01 38.37 37.07
60.13 55.09 56.78 53.86 54.36 43.50 59.04 46.47 44.08 44.66 42.33
62.68 56.61 58.28 54.51 59.89 45.00 59.55 50.89 48.25 49.93 46.78
64.47 57.09 58.72 54.15 63.96 45.30 59.10 54.36 51.53 54.15 50.31
65.93 57.05 58.63 53.24 67.46 44.95 58.19 57.36 54.61 58.05 53.55
67.12 57.36 58.91 53.07 70.40 45.33 58.00 60.14 57.49 61.71 56.57
68.36 57.37 58.83 52.60 73.02 45.18 57.24 62.59 60.04 65.06 59.27
69.04 55.76 57.33 50.35 75.21 43.55 55.03 64.77 62.26 67.93 61.67
68.86 52.54 54.39 46.09 77.09 40.40 51.58 66.64 64.02 70.32 63.65
67.63 47.83 50.20 40.41 78.48 36.24 47.11 68.18 65.59 72.23 65.31
68.51 48.30 50.54 40.65 79.76 36.83 47.38 69.69 67.27 74.01 66.96
72.96 58.97 60.33 52.79 80.76 47.68 57.62 71.06 68.79 75.60 68.33
76.48 69.94 71.10 66.09 81.63 61.33 70.30 72.18 70.03 76.82 69.47
77.78 74.79 76.09 72.42 82.13 69.00 77.14 73.09 70.96 77.80 70.43
78.27 76.83 78.20 75.23 82.23 72.99 80.45 73.88 71.78 78.69 71.28
78.30 77.76 79.19 76.59 82.10 74.93 81.86 74.38 72.29 79.31 71.90
78.20 78.01 79.39 77.00 81.80 75.81 82.38 74.87 72.80 79.85 72.56
78.38 78.23 79.61 77.23 81.86 76.28 82.60 75.37 73.28 80.53 73.14];

%% Validation
assert(size(R,1) == numel(wavelength), ...
    'R must have one row for each wavelength.');
assert(size(R,2) == numel(samples), ...
    'R must have one column for each sample.');
assert(all(R(:) > 0 & R(:) <= 100), ...
    'Reflectance values must be in the range (0, 100].');

%% Kubelka-Munk transformation
R_fraction = R / 100;
KS = (1 - R_fraction).^2 ./ (2 * R_fraction);

%% Wavelength derivative
dKS = abs(diff(KS,1,1) ./ diff(wavelength));
wavelength_d = wavelength(1:end-1);

%% Create table and save
dKS_table = array2table(dKS, 'VariableNames', samples);
dKS_table.Wavelength_nm = wavelength_d;
dKS_table = movevars(dKS_table, 'Wavelength_nm', 'Before', 1);

writetable(dKS_table, fullfile('results', 'dk_s_table.csv'));

%% Plot
figure('Color','w');
hold on;

for i = 1:numel(samples)
    plot(wavelength_d, dKS(:,i), 'LineWidth', 1.5);
end

xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|d(K/S)/d\lambda|', 'FontSize', 12, 'FontWeight', 'bold');
title('Wavelength Derivative of Kubelka-Munk K/S', ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend(samples, 'Location', 'bestoutside');
hold off;

exportgraphics(gcf, fullfile('results', 'dk_s_derivative.png'), 'Resolution', 300);

disp('K/S derivative analysis completed.');
disp(dKS_table);
