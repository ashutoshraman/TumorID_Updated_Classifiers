%Code used to normalize spectra, pegging the frequency nearest to 500 to an
%arbitrary unit of 1.0;
%Should use this code to prepare spectra for use in the simulation

load('wavelengths')
[Y, I] = min(abs((wavelengths - 500))); %I is index of wavelength nearest to 500

load('mean_healthy_spectra'); 
mean_healthy_spectra = mean_healthy_spectra .* (1. / mean_healthy_spectra(I)); %Normalize it
save('mean_healthy_spectra', 'mean_healthy_spectra'); %save it

load('mean_tumor_spectra');
mean_tumor_spectra = mean_tumor_spectra .* (1. / mean_tumor_spectra(I)); %Normalize it
save('mean_tumor_spectra', 'mean_tumor_spectra'); %save it

%plot them
%plot(wavelengths, mean_tumor_spectra, 'blue');
%hold on;
plot(wavelengths, mean_healthy_spectra, 'red');
legend({'Tumor','Healthy'}, 'FontSize', 15);
title('Spectra','FontSize',15); 
xlabel('wavelength(nm)','FontSize',15); 
ylabel('Normalized Intensity (au), Pegged 500nm at 1.0','FontSize',15);




axis([200 1000 0 1.5])
legend off
xlabel('Wavelength(nm)','FontSize',15)
ylabel('Normalized Intensity (a.u.)', 'FontSize', 15)
title('Tumor', 'FontSize', 15)
