close all;
laser_power = .175; %0.200; %Amperes, in cut_path_obj, p_points has multiplier of 80. Is this the same 80 as here
spot_size = .75; %mm  
myLaser = laser405(laser_power, spot_size); %myLaser is instance of laser405
myArgs.length = 40; %mm   % myArgs is a data struct array with length and width as args
myArgs.width = 40; %mm
myArgs.resolution = 100
mySample = sample405(myArgs); %mySample is instance of sample405 class, with myArgs struct
mySample.set_laser(myLaser); %mySample method of set_laser with myLaser instance as arg (class within class)

%set spectra
load('red_spectra'); % this is where you load red (tumor) and blue (healthy) spectra
load('blue_spectra');
load('wavelengths');
mySample.tumor_spectra = red_spectra(:, 2);
mySample.healthy_spectra = blue_spectra(:, 2);
mySample.wavelengths = wavelengths;
                                    
%Create hole for tumor and fill it w/ Tumor
length = max(size(mySample.x_coord));
width = max(size(mySample.y_coord));
thickness = mySample.thickness;

custom_lesion = create_lesion(mySample.state, length, width, thickness, 'square');
mySample.tumor_state = custom_lesion;


I = find(mySample.tumor_state == thickness);
mySample.tumor_state(:,:) = thickness;
mySample.tumor_state(I) = NaN;

step_size_array = [1, .9, .8, .7, .6, .5, .4, .3, .2, .1, .01];
int_time_array = [1, 10, 20, 30];
d_error = zeros(size(step_size_array));
for i = 1 : max(size(step_size_array))
    

    %Create a scan path
    L = 1; %length (mm)
    W = 3.6; %width (mm)
    d = 1.75; %beam width (mm)
    alpha = 1.0; %step size in terms of r between raster cuts
    dir = 'x'; %direction of cuts 'x' or 'y';
    int_time =.25; %W (mm) divided by integration time (s) gives
    step_size = step_size_array(i); %.1 step size (mm between points) is converted to points per mm through 1/step_size
    default_args = {L,W,d,alpha,dir,int_time,step_size};
    myScanPath = cut_path_obj('raster', default_args);
    loc = [10.8,20];
    myScanPath.set_location(loc); %center of mySampleopen
   
    
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
    
%     %Plot the results
%     H = figure(); %subplot(211);
%     plot_gcf(mySample, 'on', H); view(0,90);
%     hold on;
%     num_pts = max(size(myScanPath.x_points));
%     myScanPath.currentPoint = 1; %Reset count
%     k = 1;
%     while true
%         next_pt = myScanPath.nextPoint; %[t, x, y]
%         if isnan(next_pt), break; end
%         %surf(X + next_pt(2), Y + next_pt(3), beam_profile + mySample.thickness);
%         %shading interp
%         scatter3(next_pt(2), next_pt(3), max(max(beam_profile)) + mySample.thickness, 1000, 'red', 'o'); % try to make own circle
%         %scatter3(next_pt(2), next_pt(3), max(max(beam_profile)) + mySample.thickness, 10, 'red', 'o');
%         text(next_pt(2) - .3, next_pt(3) + 5, max(max(beam_profile)) + mySample.thickness, num2str(k),'FontSize',15);
%         k = k + 1;
%     end
%     axis equal
%     ylabel('Length (mm) [Y]','FontSize',15);
%     xlabel('Width (mm) [X]','FontSize',15);
    
    
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
    
   
    
    line1 = zeros(1,num_pts);
    line2 = zeros(1,num_pts);
    line3 = zeros(1,num_pts);
    x_pos = zeros(1,num_pts);
    
   
    for k = 1:num_pts
        %subplot(ceil(num_pts/2), 2,k);
        mySpectra = mySample.acquired_spectra_series(k).spectra;
        x_pos(k) = mySample.acquired_spectra_series(k).pts_xy(1);
        
        %scatter(x_pos, mySpectra(I_500) / mySpectra(I_545) - min1,'blue'); hold on; %make this more important
        line1(k) = mySpectra(I_500) / mySpectra(I_545) - min1;
        
        %scatter(x_pos, mySpectra(I_500) / mySpectra(I_575) - min2, 'red');
        line2(k) = mySpectra(I_500) / mySpectra(I_575) - min2;
        
        
        line3(k) = mySpectra(I_471) - mySpectra(I_576) - min3;
        
      
    end
    
    % Also remember to change raster graph to reflect real spot size, find
    % descending boundary using xline
    chord_length = real((2*sqrt(((length/myArgs.resolution)/4)^2 - (loc(2) - 20)^2)));
    x1_spot = 20 - chord_length/2;
    x2_spot = 20 + chord_length/2;
    
    midpt_line1 = (x_pos(find(line1==max(line1),1))+ x_pos(find(round(line1,10),1)-1))/2;
    midpt_line2 = (x_pos(find(line2==max(line2),1))+ x_pos(find(round(line2,10),1)-1))/2;
    midpt_line3 = (x_pos(find(line3==max(line3),1))+ x_pos(find(round(line3,10),1)-1))/2;
    
    d_line1 = gradient(line1,x_pos);
    d_line2 = gradient(line2,x_pos);
    d_line3 = gradient(line3,x_pos);
    
  
    predict_bounds_d_line1 = [x_pos(find(d_line1==max(d_line1))), x_pos(find(d_line1==min(d_line1)))];
    predict_bounds_d_line2 = [x_pos(find(d_line2==max(d_line2))), x_pos(find(d_line2==min(d_line2)))];
    predict_bounds_d_line3 = [x_pos(find(d_line3==max(d_line3))), x_pos(find(d_line3==min(d_line3)))];
    
  
    
    midpt_error = midpt_line1 - x1_spot;
    % d_error = predict_bounds_d_line3(1) - x1_spot
    
    d_error(i) = predict_bounds_d_line3(1) - x1_spot

end

figure();
plot(step_size_array, abs(d_error), 'k-')
xlabel('Step Size (mm between Laser Spots)')
ylabel('Predicted Boundary Error from Ground Truth (mm)')
title('Error vs. Step Size for TumorID Laser Path')

set(groot,'defaultLineLineWidth',1.2)

figure;
P1 = patch([0 x1_spot x1_spot 0],[0 0 10 10], 'cyan'); set(P1,'facealpha',0.3); hold on;
P2 = patch([x1_spot x2_spot x2_spot x1_spot],[0 0 10 10], 'magenta'); set(P2,'facealpha',0.3); hold on;

plot(x_pos, line1, 'blue'); hold on;
plot(x_pos, line2, 'green');
plot(x_pos, line3, 'red');
% plot(x_pos, line4, 'blue'); hold on;
% plot(x_pos, line5, 'red');
xline(predict_bounds_d_line1, 'yellow');

P3 = patch([x2_spot 40 40 x2_spot],[0 0 10 10], 'cyan'); set(P3,'facealpha',0.3); hold on;
%title('Tumor Edge Detection (center is tumor)','FontSize',15);
ylabel('Relative Intensity (a.u)', 'FontSize', 15); 
xlabel('X Position (mm)', 'FontSize', 15);
axis([0 40 0 0.5]);
legend({'Healthy','Tumor','500/545','500/575','471 - 576', 'Differentiation'},'FontSize',15)
axis square;
title('Step Size of .9mm')


figure();
plot(x_pos, line1); hold on;
% xline(midpt_line1)
plot(x_pos, d_line1)
xline(predict_bounds_d_line1)



return