%Code to find wavelengths with highest ratio change btwn healthy/tumor

load('mean_healthy_spectra');
load('mean_tumor_spectra');
load('wavelengths');


[Y, I_425] = min(abs((mySample.wavelengths - 425))); %I_x is index of wavelength nearest to x
[Y, I_700] = min(abs((mySample.wavelengths - 700))); %I_x is index of wavelength nearest to x

mean_healthy_spectra = mean_healthy_spectra(I_425:I_700);
mean_tumor_spectra = mean_tumor_spectra(I_425:I_700);
wavelengths = wavelengths(I_425:I_700);

num_pts = numel(wavelengths)

max_ratio = -1;
max_Inum = -1;
max_Idenom = -1;

for k = 1:num_pts
    ratios_healthy = mean_healthy_spectra(k) ./ mean_healthy_spectra;
    ratios_tumor = mean_tumor_spectra(k) ./ mean_tumor_spectra;
    
    
    %diff_ratios = abs((ratios_healthy - ratios_tumor)) ./ max([ratios_tumor; ratios_healthy]);
    diff_ratios = (mean_healthy_spectra(k) - mean_healthy_spectra) - (mean_tumor_spectra(k) - mean_tumor_spectra);
    
    %diff_ratios = abs(mean_healthy_spectra(k) - mean_healthy_spectra) - abs(mean_tumor_spectra(k) - mean_tumor_spectra);
    
    [diff_ratio Idenom] = max(abs(diff_ratios));
    
    %save highest
    if diff_ratio > max_ratio
        max_ratio = diff_ratio;
        max_Inum = k;
        max_Idenom = Idenom;
    end      
end

max_ratio
max_Inum
max_Idenom

wave_num = wavelengths(max_Inum)
wave_den = wavelengths(max_Idenom)

mean_healthy_spectra(max_Inum) / mean_healthy_spectra(max_Idenom) - mean_tumor_spectra(max_Inum) / mean_tumor_spectra(max_Idenom)