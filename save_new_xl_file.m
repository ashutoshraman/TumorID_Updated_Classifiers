blue_spectra = xlsread('blue.csv.csv');
% save red_spectra intensity
figure;
plot(blue_spectra(:, 1), blue_spectra(:, 2))
save blue_spectra.mat blue_spectra
