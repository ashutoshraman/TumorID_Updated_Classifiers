%%Code for Figure showing Tumor vs Healthy Spectra

close all;
%set spectra
load('mean_healthy_spectra');
load('mean_tumor_spectra');
load('wavelengths');

figure; 
plot(wavelengths, smoothed_healthy_spectra(1:end-7)); hold on;
plot(wavelengths, smoothed_tumor_spectra(1:end-7));
plot(wavelengths, smoothed_tumor_spectra(1:end-7) - smoothed_healthy_spectra(1:end-7) + 1.5);
line([0,1000],[1.5,1.5],'color','black');
legend({'Healthy','Tumor','(Tumor - Healthy)'},'FontSize', 15);
xlabel('Wavelengths','FontSize', 15);
ylabel('Intensities (normalized)','FontSize', 15);
legend boxoff
axis([200 1000 0 1.7])

%%


