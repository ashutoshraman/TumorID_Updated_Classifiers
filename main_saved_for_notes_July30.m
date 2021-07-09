%Code to simulate Matt's laser moving across the surface of a sample
close all;
laser_power = 80; %0.200; %Amperes
spot_size = 3.0; %mm
myLaser = laser405(laser_power, spot_size);
myArgs.length = 40; %mm
myArgs.width = 40; %mm
mySample = sample405(myArgs);
mySample.set_laser(myLaser);

%set spectra
load('mean_healthy_spectra');
load('mean_tumor_spectra');
load('wavelengths');
mySample.tumor_spectra = mean_tumor_spectra;
mySample.healthy_spectra = mean_healthy_spectra;
mySample.wavelengths = wavelengths;
                                    
%Create hole for tumor and fill it w/ Tumor
length = max(size(mySample.x_coord));
width = max(size(mySample.y_coord));
thickness = mySample.thickness;
mySample.state(length/4:length/4 + length/2, width/4:width/4 + width/2) = thickness / 2; %carve out square of half of depth
mySample.tumor_state = mySample.state;
I = find(mySample.tumor_state == thickness);
mySample.tumor_state(:,:) = thickness;
mySample.tumor_state(I) = NaN;

%Create a scan path
L = 1; %length (mm)
W = 30; %width (mm)
d = 1.75; %beam width (mm)
alpha = 1.0; %step size in terms of r between raster cuts
dir = 'x'; %direction of cuts 'x' or 'y';
speed = 1; %speed of cut mm/s
pp_cm = 1; %points per centimeter (actually mm, maybe fix this notation)
default_args = {L,W,d,alpha,dir,speed,pp_cm};
myScanPath = cut_path_obj('raster', default_args);
myScanPath.set_location([20,20]); %center of mySample

% %Create a scan path
% myScanPath = cut_path_obj('empty', -1);
% myScanPath.set_x_points([5  10 20 20 28 32]);
% myScanPath.set_y_points([20 10 20 30 25 15]);
% myScanPath.p_points = [80 80 80 80 80 80];
% myScanPath.t_points = [0.1 0.2 0.3 0.4 0.5 0.6];
% %myScanPath.set_location([0,0]); %center of mySample

%Pre allocate space for acquired spectra
num_pts = max(size(myScanPath.t_points));
sz_wavelength = max(size(mySample.wavelengths));
mySample.acquired_spectra_series = repmat(struct('pts_xy',[NaN, NaN], 'spectra', nan(1,sz_wavelength)), num_pts, 1);                                    

%execute the path
mySample.perform_ablation(myScanPath); %not actually ablating, just same function which calls acquire spectra function

%Generate a laser spot
lspace = linspace(-1.5*myLaser.spot_size,1.5*myLaser.spot_size,99);
[X,Y] = meshgrid(lspace(1:2:end),lspace(1:2:end));
beam_profile = myLaser.pwr_function(sqrt((X.^2 + Y.^2)));
beam_profile(beam_profile < 0.01) = NaN;
beam_profile = beam_profile ./ max(max(beam_profile)); %normalize for plotting
%surf(beam_profile); axis square;

%Plot the results
H = subplot(211);
plot_gcf(mySample, 'on', H); view(0,90);
hold on;
num_pts = max(size(myScanPath.x_points));
myScanPath.currentPoint = 1; %Reset count
while true
    next_pt = myScanPath.nextPoint; %[t, x, y]
    if isnan(next_pt), break; end
    surf(X + next_pt(2), Y + next_pt(3), beam_profile + mySample.thickness);
    shading interp
end
axis equal

%Smooth healthy and tumor spectra for comparison
smoothed_healthy_spectra = conv(mySample.healthy_spectra, [1 1 1 1 1 1 1 1] ./ 8);
smoothed_tumor_spectra = conv(mySample.tumor_spectra, [1 1 1 1 1 1 1 1] ./ 8);

H = figure; %for the spectra
for k = 1:num_pts
    %subplot(ceil(num_pts/2), 2,k);
    figure(H); H.clo;
    plot(mySample.wavelengths, smoothed_healthy_spectra(8:end), 'k');
    hold on;
    plot(mySample.wavelengths, smoothed_tumor_spectra(8:end), 'g');
    smoothed_spectra = conv(mySample.acquired_spectra_series(k).spectra, [1 1 1 1 1 1 1 1] ./ 8);
    plot(mySample.wavelengths, smoothed_spectra(8:end), '--r');
    title(sprintf('Spectra %d', k));
    ylabel('Normalized Itensity A.U. (500nm -> 1.0)', 'FontSize', 15); 
    legend({'Healthy','Tumor','Acquired'},'FontSize',15);    
    xlabel('Wavelengths (nm)', 'FontSize', 15);
    axis([0 1200 0 1.5]);
    pause
end

[Y, I_500] = min(abs((mySample.wavelengths - 500))); %I_x is index of wavelength nearest to x
[Y, I_545] = min(abs((mySample.wavelengths - 545))); %I_x is index of wavelength nearest to x
[Y, I_575] = min(abs((mySample.wavelengths - 575))); %I_x is index of wavelength nearest to x
[Y, I_425] = min(abs((mySample.wavelengths - 425))); %I_x is index of wavelength nearest to x
[Y, I_750] = min(abs((mySample.wavelengths - 750))); %I_x is index of wavelength nearest to x

H = gcf;
figure(H);
subplot(212); %for the spectra
for k = 1:num_pts
    %subplot(ceil(num_pts/2), 2,k);
    mySpectra = mySample.acquired_spectra_series(k).spectra;
    x_pos = mySample.acquired_spectra_series(k).pts_xy(1);
    scatter(x_pos, mySpectra(I_500) / mySpectra(I_545),'blue'); hold on; %make this more important
    scatter(x_pos, mySpectra(I_500) / mySpectra(I_575), 'red');
    scatter(x_pos, mean(mySpectra(I_425:I_750)),'green');
    hold on;
    title('Tumor Edge Detection (center is tumor)','FontSize',15);
    ylabel('Relative Intensity', 'FontSize', 15); 
    xlabel('X Position (mm)', 'FontSize', 15);
    axis([0 40 0 2.5]);
    legend({'500/545','500/575','Mean Intensity 425~750'},'FontSize',15)
    axis square;
end

return
