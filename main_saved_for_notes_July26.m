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

% %Create a scan path
% L = 1; %length (mm)
% W = 30; %width (mm)
% d = 1.75; %beam width (mm)
% alpha = 1.0; %step size in terms of r between raster cuts
% dir = 'x'; %direction of cuts 'x' or 'y';
% speed = 10; %speed of cut mm/s
% pp_cm = 0.2; %points per centimeter (actually mm, maybe fix this notation)
% default_args = {L,W,d,alpha,dir,speed,pp_cm};
% myScanPath = cut_path_obj('raster', default_args);
% myScanPath.set_location([22,20]); %center of mySample

%Create a scan path
myScanPath = cut_path_obj('empty', -1);
myScanPath.set_x_points([5  10 20 20 28 32]);
myScanPath.set_y_points([20 10 20 30 25 15]);
myScanPath.p_points = [80 80 80 80 80 80];
myScanPath.t_points = [0.1 0.2 0.3 0.4 0.5 0.6];
%myScanPath.set_location([0,0]); %center of mySample

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
plot(mySample); view(0,90);
hold on;
num_pts = max(size(myScanPath.x_points));
myScanPath.currentPoint = 1; %Reset count
k = 1;
while true
    next_pt = myScanPath.nextPoint; %[t, x, y]
    if isnan(next_pt), break; end
    surf(X + next_pt(2), Y + next_pt(3), beam_profile + mySample.thickness);
    text(next_pt(2) - 0.5, next_pt(3), max(max(beam_profile)) + mySample.thickness, num2str(k),'FontSize',15,'Color','red');
    k = k + 1;
end
axis equal

%Smooth healthy and tumor spectra for comparison
smoothed_healthy_spectra = conv(mySample.healthy_spectra, [1 1 1 1 1 1 1 1] ./ 8);
smoothed_tumor_spectra = conv(mySample.tumor_spectra, [1 1 1 1 1 1 1 1] ./ 8);

H = figure; %for the spectra
subplot(321);
for k = 1:num_pts
    %subplot(ceil(num_pts/2), 2,k);
    %figure(H); H.clo;
    subplot(2,3,k);
    plot(mySample.wavelengths, smoothed_healthy_spectra(8:end), 'k');
    hold on;
    plot(mySample.wavelengths, smoothed_tumor_spectra(8:end), 'g');
    smoothed_spectra = conv(mySample.acquired_spectra_series(k).spectra, [1 1 1 1 1 1 1 1] ./ 8);
    plot(mySample.wavelengths, smoothed_spectra(8:end), '--r');
    %title(sprintf('Spectra %d', k));
    ylabel('Normalized Intensity', 'FontSize', 15); 
    if (k == 3)
        legend({'Healthy','Tumor','Acquired'},'FontSize',13); legend boxoff;
    end
    xlabel('Wavelengths (nm)', 'FontSize', 15);
    text(300,1.2,join(['(',num2str(k),')']),'FontSize',15);
    axis([200 1200 0 1.5]);
end
