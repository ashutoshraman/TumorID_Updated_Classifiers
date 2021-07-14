red_spectra = xlsread('red.csv.csv');
% save red_spectra intensity
[Y, I_500] = min(abs((mySample.wavelengths - 500))); %I_x is index of wavelength nearest to x
red_spectra = [red_spectra(:, 1), red_spectra(:, 2)./red_spectra(I_500, 2)]
figure;
plot(red_spectra(:, 1), red_spectra(:, 2))
% plot(red_spectra(:, 1), red_spectra(:, 2)./red_spectra(I_500, 2))
save red_spectra.mat red_spectra
