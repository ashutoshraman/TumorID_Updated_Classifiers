%Code to simulate Matt's laser moving across the surface of a sample
%Intent is to image the entire sample and then generate an image of the
%sample from the recorded spectra only.

close all; clearvars;

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
% L = 30; %length (mm)
% W = 30; %width (mm)
% d = 0.5; %beam width (mm)
% alpha = 1.0; %step size in terms of r between raster cuts
% dir = 'x'; %direction of cuts 'x' or 'y';
% speed = 1; %speed of cut mm/s
% pp_mm = 0.5; %points per millimeter
% default_args = {L,W,d,alpha,dir,speed,pp_mm};
% %myScanPath = cut_path_obj('raster', default_args);
% %myScanPath.set_location([20,20]); %center of mySampleopen

%Create a scan path
myScanPath = cut_path_obj('empty', -1);
pts = -15:2:15;
[X,Y] = meshgrid(pts,pts);
myScanPath.set_x_points(reshape(X,1,numel(X)));
myScanPath.set_y_points(reshape(Y,1,numel(Y)));
myScanPath.p_points = 80*ones(size(myScanPath.x_points));
myScanPath.t_points = 0.1:0.1:0.1*numel(myScanPath.x_points);
myScanPath.set_location([20,20]); %center of mySampleopen

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


%make a circle shape for plotting
r_spot = spot_size * mySample.resolution / 2;
theta = 0:0.1:2.1*pi;
x = cos(theta)*spot_size/2;
y = sin(theta)*spot_size/2;

%Plot the results
H = figure(); %subplot(211);
plot_gcf(mySample, 'on', H); view(0,90);
hold on;
num_pts = max(size(myScanPath.x_points));
myScanPath.currentPoint = 1; %Reset count
k = 1;
while true
    next_pt = myScanPath.nextPoint; %[t, x, y]
    if isnan(next_pt), break; end
    plot3(x+next_pt(2), y+next_pt(3), mySample.thickness*ones(size(x)),'green');
    k = k + 1;
end
axis equal
ylabel('Length (mm) [Y]','FontSize',15);
xlabel('Width (mm) [X]','FontSize',15);
title('Acquisition Spots','FontSize',15);

drawnow

%Smooth healthy and tumor spectra for comparison
smoothed_healthy_spectra = conv(mySample.healthy_spectra, [1 1 1 1 1 1 1 1] ./ 8);
smoothed_tumor_spectra = conv(mySample.tumor_spectra, [1 1 1 1 1 1 1 1] ./ 8);

[Y, I_500] = min(abs((mySample.wavelengths - 500))); %I_x is index of wavelength nearest to x
[Y, I_545] = min(abs((mySample.wavelengths - 545))); %I_x is index of wavelength nearest to x
[Y, I_575] = min(abs((mySample.wavelengths - 575))); %I_x is index of wavelength nearest to x

%max ratio change between healthy and tumor between 475 and 700nm only
[Y, I_527] = min(abs((mySample.wavelengths - 527))); %I_x is index of wavelength nearest to x
[Y, I_687] = min(abs((mySample.wavelengths - 687))); %I_x is index of wavelength nearest to x

%max intensity change between healthy and tumor between 475 and 700nm only
[Y, I_527] = min(abs((mySample.wavelengths - 527))); %I_x is index of wavelength nearest to x
[Y, I_576] = min(abs((mySample.wavelengths - 576))); %I_x is index of wavelength nearest to x

%max difference 
[Y, I_471] = min(abs((mySample.wavelengths - 471))); %I_x is index of wavelength nearest to x
[Y, I_576] = min(abs((mySample.wavelengths - 576))); %I_x is index of wavelength nearest to x


mySpectra = mySample.acquired_spectra_series(1).spectra;

min1 = mySpectra(I_500) / mySpectra(I_545);
min2 = mySpectra(I_500) / mySpectra(I_575);
min3 = mySpectra(I_471) - mySpectra(I_576);

%Make empty "image" for the results to go into
image_reconstruction = zeros(size(mySample.state));
image_reconstruction_sample_number = zeros(size(image_reconstruction));
%h_reconstruction = figure();

myScanPath.currentPoint = 1; %Reset count
k = 1;
while true
    next_pt = myScanPath.nextPoint; %[t, x, y]
    if isnan(next_pt), break; end
    mySpectra = mySample.acquired_spectra_series(k).spectra;
    x_pos = next_pt(2);
    y_pos = next_pt(3);
    
    classifier_value = mySpectra(I_500) / mySpectra(I_545);
    
    x_temp = mySample.resolution * x_pos;
    y_temp = mySample.resolution * y_pos;
    r_spot = spot_size * mySample.resolution / 2;
    
    rad_steps = ceil(r_spot);
    
    x_coords_near = x_temp - rad_steps : x_temp + rad_steps;
    y_coords_near = y_temp - rad_steps : y_temp + rad_steps;
    [X,Y] = meshgrid(x_coords_near, y_coords_near);
    
    usable_pts = sqrt((X - x_temp).^2 + (Y - y_temp).^2) < rad_steps;
    %usable_pts = sqrt((x_coords_near - x_temp).^2 + (y_coords_near - y_temp).^2) < rad_steps;
    
    x_coords_near = X(usable_pts);
    y_coords_near = Y(usable_pts);
    
    index_near = sub2ind(size(image_reconstruction_sample_number),x_coords_near, y_coords_near);
    
    %Make reconstruction the average of each overlapping area
    image_reconstruction_sample_number(index_near) = image_reconstruction_sample_number(index_near) + 1; %track how many samples are taken of each spot.
    scaling_ratio = 1 ./ image_reconstruction_sample_number(index_near);
    scaling_ratio(scaling_ratio == inf) = 0;
    image_reconstruction(index_near) = scaling_ratio.*classifier_value + (1-scaling_ratio).* image_reconstruction(index_near);
    
    if classifier_value > .90
        %figure(h_reconstruction); hold on;
        %scatter(x_pos, y_pos, 2000, 'red', 'o');
    else
        %figure(h_reconstruction); hold on;
        %scatter(x_pos, y_pos, 2000, 'black', 'o');
    end
    k = k + 1;
end
axis([0 40 0 40]);

drawnow;


h1 = figure()
imshow(image_reconstruction ./ max(max(image_reconstruction)));
colormap(gca, jet(256)); caxis([0.8 1.0]);
rectangle('Position', [mySample.length/4 mySample.width/4 mySample.length/2 mySample.width/2]*mySample.resolution,'EdgeColor', 'Black', 'LineWidth', 2, 'LineStyle', '--');
title(h1.CurrentAxes, 'Reconstruction','FontSize',15)
axis equal;
ylabel('Length (mm) [Y]','FontSize',15);
xlabel('Width (mm) [X]','FontSize',15);


% h2 = figure()
% imshow(image_reconstruction_sample_number ./ max(max(image_reconstruction_sample_number))); 
% title(h2.CurrentAxes, 'num samples')
return