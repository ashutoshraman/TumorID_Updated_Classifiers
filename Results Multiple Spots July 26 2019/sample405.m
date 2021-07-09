classdef sample405 < matlab.mixin.Copyable
    properties
        thickness
        width %x direction
        length %y direction
        resolution
        
        state %surface created from max(tumor_state,healthy_state)
        tumor_state %surface of tumor
        healthy_state %surface of healthy tissue
        
        histology
        x_coord %coordinates in x (width) direction
        y_coord %coordinates in y (length) direction
        laser %for ablatoin
        laser405 %for interrogation
        save_images = false;
        plot_tissue = false;
        folder_name %folder where frames are saved for video making
        %these two are used for fixed parametes from A4 fiting
        alpha %density * h_abl
        irrad_thresh %for fixed irradiance thresholds
        
        total_irradiance %for Matt's stuff. Sum of total irradiance, irreguardless of irrad_thresh
        %material ablation properties
        %Theta_th %Threshold radiant exposure for material ejection
        %density %density g/mm^3
        %h_abl %ablation enthalpy J/g
        %mu_a %absorption coefficient, mm^-1
        material_properties %cell of nxm tissue properties (density, h_abl, and mu_a)
                %theta_th_xy %calculated thetas for each position on the surface of the tissue.
                %density %density g/mm^3
                %h_abl %ablation enthalpy J/g
                %mu_a %absorption coefficient, mm^-1
        r_th %radius of ablation zone for the given laser
        
        %Spectra
        tumor_spectra %normalized
        healthy_spectra %normalized
        wavelengths
        
        %Data out
        acquired_spectra_series %matrix of acquired spectra, for each point in the blah blah
    end
    
    methods
        function obj = sample405(varargin)
            myArgs = varargin{1};
            if isfield(myArgs,'width')
                obj.width = myArgs.width;
            else
                obj.width = 10;
                disp('Sample405: Setting sample width to 10mm (default).');
            end %
            
            if isfield(myArgs,'length')
                obj.length = myArgs.length;
            else
                obj.length = obj.width;
                disp('Sample405: Setting sample length to sample width (default).');
            end %
            
            if isfield(myArgs,'alpha')
                obj.alpha = myArgs.alpha;
            else
                obj.alpha = 3.3285;
                disp('Sample405: Setting sample alpha to 3.3285 (default).');
            end %
            
            if isfield(myArgs,'irrad_thresh')
                obj.irrad_thresh = myArgs.irrad_thresh;
            else
                obj.irrad_thresh = 0.3107;
                disp('Sample405: Setting sample irrad_thresh to 0.3107 (default).');
            end %
           
            if isfield(myArgs,'thickness')
                obj.thickness = myArgs.thickness;
            else
                obj.thickness = 5; %mm
                disp('Sample405: Setting sample thickness to 5mm (default).');
            end %
            
            if isfield(myArgs,'resolution')
                obj.resolution = myArgs.resolution;
            else
                obj.resolution = 10; %number of points per mm
                disp('Sample405: Setting sample resolution to 10 (default).');
            end
            
            
            if all(isfield(myArgs,{'state','x_coord','y_coord'}))
                obj.state = myArgs.state;
                obj.x_coord = myArgs.x_coord;
                obj.y_coord = myArgs.y_coord;
            else
                obj.x_coord = linspace(0, obj.width, obj.width*obj.resolution);
                obj.y_coord = linspace(0, obj.length, obj.length*obj.resolution);
                obj.state = obj.thickness*ones(obj.width*obj.resolution, obj.length*obj.resolution);
                disp('Sample405: Setting state (default).');
            end
            
            if isfield(myArgs,'tumor_state')
                obj.tumor_state = myArgs.tumor_state;
            else
                obj.tumor_state = zeros(size(obj.state)); %number of points per mm
                disp('Sample405: Setting tumor state to zero plane (default).');
            end
            
            %This is for Matt's stuff. 
            obj.total_irradiance = zeros(size(obj.state));
            
            obj.alpha = obj.alpha*ones(size(obj.state));  %from fitting on A4
                %split the sample in half by properties
                %[La Wa] = size(obj.alpha); 
                %obj.alpha(La/2:end, :) = obj.alpha(La/2:end,:) * 2;
            
            obj.irrad_thresh = obj.irrad_thresh*ones(size(obj.state));
            %[a,b] = size(obj.state);
            %obj.irrad_thresh(:,b/2:end) = obj.irrad_thresh(a:end, b:end) .* 20;
            
            %obj.histology = load('tumor_tissue');
            %obj.density = 1 / 1000; %g/mm^3  (1.073g/cm^3 for human frontal white, 1.090 for human frontal gray, "The density of tissues in and about the head", 1970, Barber
            %obj.h_abl = 2580; %J/g
            %obj.mu_a = 100; %cm^-1
            %obj.Theta_th = obj.density * obj.h_abl / obj.mu_a;
            
            %Material Properties for each voxel, rather than homogenous
            %density_white = 1.073/1000; %g/mm^3
            %density_gray = 1.073/1000; %g/mm^3
            
            %obj.material_properties.density = density_white;
%             obj.material_properties.density = ones(obj.width*obj.resolution, obj.length*obj.resolution); 
%                 [L W] = size(obj.material_properties.density);
%                 obj.material_properties.density(1:L/2,:) = density_white;
%                 obj.material_properties.density(L/2:end,:) = density_white;
            
            
            %obj.material_properties.h_abl = 2000;
            %obj.material_properties.mu_a = obj.irrad_thresh / obj.material_properties.density / obj.material_properties.h_abl;
            %obj.material_properties.theta_th_xy = obj.material_properties.density * obj.material_properties.h_abl * obj.material_properties.mu_a;      
            
%             The problem right now is that with the new laser fitting, the absolute power is much less than what I've been testing with. 
%             now the incident irradiation does not ever get above the threshold. Fix this.
            
            %obj.laser = laser();            
            %                 %laser, replace later with a laser obj
            %                 obj.laser.width = 1;
            %                 obj.laser.position.x = 2; %this is center of beam in x (mm)
            %                 obj.laser.position.y = 2; %this is center of beam in y (mm)
            %                 obj.laser.cutting_rate = 1; %dh/dt
            %                 obj.laser.power = 'on'; %or 'off'
            %                 lspace = linspace(-obj.laser.width/2,obj.laser.width/2,21);
            %                 [X,Y] = meshgrid(lspace,lspace);
            %                 obj.laser.beam_profile = exp(-10*(X.^2 + Y.^2))
        end
        
        function update_ablation_steady_state_model(obj, dt)
            %this function updates the ablation using material parameters
            %and irradiance values.  This is the more accurate one.
            %We will calculate the dz for this time step at every
            %location on the surface of the sample
            
            %find indicies of location of laser on surface
            x_centers = obj.x_coord + 0.5/obj.resolution;
            x_centers = x_centers(x_centers < obj.width);
            
            y_centers = obj.y_coord + 0.5/obj.resolution;
            y_centers = y_centers(y_centers < obj.width);
            
            [xdist, Ix] = min(abs(x_centers-obj.laser.position.x));
            [ydist, Iy] = min(abs(y_centers-obj.laser.position.y));
            
            
            vox_area = (1/obj.resolution)^2;

            rad_steps = ceil(obj.r_th * obj.resolution);
            
            x_coords_near = Ix - rad_steps : Ix + rad_steps;
            y_coords_near = Iy - rad_steps : Iy + rad_steps;
            usable_pts = (x_coords_near > 0) & (y_coords_near > 0) & (x_coords_near < obj.width * obj.resolution) & (y_coords_near < obj.length * obj.resolution);
            x_coords_near = x_coords_near(usable_pts);
            y_coords_near = y_coords_near(usable_pts);
            
            irrad_fun = obj.laser.pwr_function;
           
            Theta_th_near = obj.irrad_thresh(x_coords_near, y_coords_near); %J/mm^2
            alpha_near = obj.alpha(x_coords_near, y_coords_near); %J/mm^3
            
            centers_x = x_centers(x_coords_near);
            x_dists = centers_x - obj.laser.position.x;
            
            centers_y = y_centers(y_coords_near);
            y_dists = centers_y - obj.laser.position.y;
            
            [X,Y] = meshgrid(x_dists, y_dists); %Flipped this, maybe fixes our rotation issue?
            irrad_levels = irrad_fun(sqrt(X.^2 + Y.^2)); %W/mm^2
            incident_irradiation = dt.*irrad_levels; %For Matt's stuff
            obj.total_irradiance(x_coords_near, y_coords_near) = obj.total_irradiance(x_coords_near, y_coords_near) + incident_irradiation; %for Matt's stuff
            
            
            Thetas_0 = irrad_levels * 1; %W/mm^2; assume a 1s pulse, then scale depth by time.
            dzs = dt*max(0,(Thetas_0 - Theta_th_near)) ./ alpha_near;
            
            %Original Code
            %obj.state(x_coords_near, y_coords_near) = obj.state(x_coords_near, y_coords_near) - dzs;
            
            %Code to work with tumor and healthy tissue
            if ~isempty(obj.tumor_state)
                mask = (obj.tumor_state >= obj.state);
                obj.tumor_state(~mask) = nan;
                obj.tumor_state(x_coords_near, y_coords_near) = obj.tumor_state(x_coords_near, y_coords_near) - mask(x_coords_near, y_coords_near).*dzs;
                obj.state(x_coords_near, y_coords_near) = obj.state(x_coords_near, y_coords_near) - ~mask(x_coords_near, y_coords_near).*dzs;
            else
                obj.state(x_coords_near, y_coords_near) = obj.state(x_coords_near, y_coords_near) - dzs;
            end
        end       
        
        function acquired_spectra = acquire_spectra(obj, dt)
            %function to acquire spectral signature based on energy in
            %and known spectral signature for a unit energy input
            %output is a spectral signature, which is an integration over
            %the spot size of the interrogation laser=
            
            %find indicies of location of laser on surface
            x_centers = obj.x_coord + 0.5/obj.resolution;
            x_centers = x_centers(x_centers < obj.width);
            
            y_centers = obj.y_coord + 0.5/obj.resolution;
            y_centers = y_centers(y_centers < obj.width);
            
            [xdist, Ix] = min(abs(x_centers-obj.laser.position.x));
            [ydist, Iy] = min(abs(y_centers-obj.laser.position.y));
            
            vox_area = (1/obj.resolution)^2;

            rad_steps = ceil(obj.r_th * obj.resolution);
            
            x_coords_near = Ix - rad_steps : Ix + rad_steps;
            y_coords_near = Iy - rad_steps : Iy + rad_steps;
            usable_pts = (x_coords_near > 0) & (y_coords_near > 0) & (x_coords_near < obj.width * obj.resolution) & (y_coords_near < obj.length * obj.resolution);
            x_coords_near = x_coords_near(usable_pts);
            y_coords_near = y_coords_near(usable_pts);
            
            irrad_fun = obj.laser.pwr_function;
            
            centers_x = x_centers(x_coords_near);
            x_dists = centers_x - obj.laser.position.x;
            
            centers_y = y_centers(y_coords_near);
            y_dists = centers_y - obj.laser.position.y;
            
            [X,Y] = meshgrid(x_dists, y_dists); %Flipped this, maybe fixes our rotation issue?
            irrad_levels = irrad_fun(sqrt(X.^2 + Y.^2)); %W/mm^2
            incident_irradiation = dt.*irrad_levels; %For Matt's stuff
            obj.total_irradiance(x_coords_near, y_coords_near) = obj.total_irradiance(x_coords_near, y_coords_near) + incident_irradiation; %for Matt's stuff 
            
            %# out @ wavelength = intensity * area * seconds * spectral_intensity(wavelength)
            %Spectra Out = sum(intensity * area * seconds *
            %specral_intensity(wavelength)) for all wavelengths
            %Spectra Out = area*seconds*
            %n x n * n * n * num(wavelength) = x * num(wavelength)
            
            %Generate labeled (tumor/healthy) surface 
            near_healthy = obj.state(x_coords_near, y_coords_near);
            near_tumor = obj.tumor_state(x_coords_near, y_coords_near);
            I_t = find(near_tumor >= near_tumor); %I_t are the indices of the tumor
            labeled_surface = zeros(size(near_healthy));
            labeled_surface(I_t) = 1; %Label tumor as 1's
            
            %Calculate energy / photons in per tissue type
            energy_in_tumor = sum(sum(labeled_surface .* incident_irradiation));
            energy_in_healthy = sum(sum(~labeled_surface .* incident_irradiation));
            
            %calculated "sensed" spectra as sum of contributions from
            %healthy and tumor tissue
            acquired_spectra = energy_in_tumor * obj.tumor_spectra + energy_in_healthy * obj.healthy_spectra;
            
            %normalize
            [Y, I] = min(abs((obj.wavelengths - 500))); %I is index of wavelength nearest to 500
            acquired_spectra = transpose(acquired_spectra .* (1. / acquired_spectra(I))); %Normalize and transpose it
        end
        
        function perform_ablation(obj, cut_path)%material,laser,dt)
            %path = [t; laser_x; laser_y], size = 3 x N
            %dt = 0.01; %cut_path(2,1) - cut_path(1,1);   
%             t = cut_path(1,:);
%             dt = diff(t); %time BETWEEN points. size = n - 1;
%             dt = [dt, dt(end)]; %Repeat last dt for last point. I don't have a better method.
            obj.r_th = obj.calculate_max_ablation_r();
            
            t_final = cut_path.t_points(end);
            t_old = 0;
            t_new = 0;
            t_last_frame = -1.0;
            t_last_printout = -1.0;
            
            %[L W] = size(cut_path.t);
            for k = 1:length(cut_path.t_points)
                %dt = cut_path(1,k);
%                 obj.laser.position.x = cut_path(2,k);
%                 obj.laser.position.y = cut_path(3,k);
%                 obj.update_ablation_steady_state_model(dt(k));
                
                %nextPoint = cut_path.nextPoint();
                nextPoint = cut_path.getPoint(k);
                
                t_new = nextPoint(1);
                obj.laser.position.x = nextPoint(2);
                obj.laser.position.y = nextPoint(3);
                
                %figure(5); hold on; scatter(nextPoint(2), nextPoint(3));
                %figure(4);
                
                %obj.update_ablation_steady_state_model(t_new - t_old);
                acquired_spectra = obj.acquire_spectra(t_new - t_old); %Added for spectral acquisition
                
                %Save the spectra and x/y pt of laser to acquired_spectra
                obj.acquired_spectra_series(k).pts_xy = [obj.laser.position.x, obj.laser.position.x];
                obj.acquired_spectra_series(k).spectra = acquired_spectra;

                
                %temporary code
                %H_spectra = figure(); plot(obj.wavelengths, acquired_spectra); title('Acquired');
                %pause
                %close(H_spectra);
                
                %save images for video
                if obj.save_images && (t_new - t_last_frame > 1/30) %capture at frame rate of 30
                    t_last_frame = t_new; 
                    obj.plot_gcf('off');
                    %img = getframe(gcf);
                    %imwrite(img.cdata,fullfile('./',obj.folder_name,sprintf('img%d.jpg',k)));
                    print(gcf,fullfile('./',obj.folder_name,sprintf('img%d.jpg',k)),'-dpng','-r300');
                end
                
                
                if ((mod(k,10) == 1) && obj.plot_tissue) %plot only every 10 steps, for speed
                    obj.plot_gcf('on');
                    pause(0.1);
                end
                
                %keep taps on simulation length
                if t_new - t_last_printout > 1
                    t_last_printout = t_new;
                    %disp(sprintf('t=%f of %f (simulated time)',t_new,t_final));
                end
                
                t_old = t_new; %update time.

            end
        end
        
        function set_save_images(obj, value)
            obj.save_images = value; %should be boolean
            obj.folder_name = date();
            if ~isdir(obj.folder_name)
                mkdir(obj.folder_name);
            end
        end
        
        function set_plot_tissue(obj, value)
            obj.plot_tissue = value; %should be boolean
        end
        
        function xs = get_cross_section_x(obj, y_coord)
            %takes cross section in the X direction
            %at a specified Y coordinate.
            xs = obj.state(y_coord*obj.resolution, :);
        end
        
        function xs = get_cross_section_y(obj, x_coord)
            %takes cross section in the Y direction
            %at a specified X coordinate.
            xs = obj.state(:, x_coord*obj.resolution)';
        end
        
        function state = get_tissue_state(obj, xy_range)
            %alpha = obj.resolution;
            xmin = xy_range(1)*obj.resolution;
            xmax = xy_range(2)*obj.resolution;
            ymin = xy_range(3)*obj.resolution;
            ymax = xy_range(4)*obj.resolution;
            
            state = obj.state(xmin:xmax, ymin:ymax);
        end
        
        function set_laser(obj, laser)
            obj.laser = laser;
        end
        
        function plot(obj)
            H = figure;
            obj.plot_gcf('on', H)
        end
        
        function plot_gcf(obj,visibility, fig_handle)
            
            %%%Original Code
            if (false)
                h = gcf();
                set(h,'visible',visibility);
                [az,el] = view;  %Save current perspective
                h.clo; %clear

                [X,Y] = meshgrid(obj.x_coord, obj.y_coord);
                h = surf(X, Y, obj.state);

                view(az,el);
                set(h,'LineStyle','none')
                axis equal;
                axis([0 obj.width 0 obj.length 0 obj.thickness*1.5]);
                %plot laser
                %scatter3(obj.laser.position.y, obj.laser.position.x, obj.thickness*1.2, 50, 'O', 'MarkerEdgeColor', 'black');
                %axis equal;  
                xlabel('Width [mm] (X)','FontSize',13); 
                ylabel('Length [mm] (Y)','FontSize',13);

                caxis([floor(min(obj.state(:))), ceil(max(obj.state(:)))]);
            end
            %%End Original Code
            
            
            %%%%%%%%
            
            %%New code plotting tumor and healthy
            h = fig_handle;
            set(h,'visible',visibility);
            [az,el] = view;  %Save current perspective
            h.clo; %clear

            [X,Y] = meshgrid(obj.x_coord, obj.y_coord);
            mask = (obj.tumor_state >= obj.state);
            mySurf = max(obj.tumor_state, obj.state);
            h = surf(X, Y, mySurf, mask+3);
            
            view(az,el);
            set(h,'LineStyle','none')
            axis equal;
            axis([0 obj.width 0 obj.length 0 obj.thickness*1.5]);
            
            
            %plot laser
            %Todo: Make this representative of the region of effect of the
            %laser
            % hold on; scatter3(obj.laser.position.y, obj.laser.position.x, obj.thickness*1.2, 50, 'O', 'MarkerEdgeColor', 'black');

            xlabel('Width [mm] (X)','FontSize',13); 
            ylabel('Length [mm] (Y)','FontSize',13);
            
            %caxis([floor(min(obj.state(:))), ceil(max(obj.state(:)))]);
            
            %cmap = colormap();
            %colormap([255 255 183]./255); %like skin
            light; %add light for shading

            %%end new code plotting tumor and healthy
            
            %%%%%%%%
            
            %%%%%%%Plot tissue sides
            %cmap = colormap();
            %colormap([255 255 183]./255); %like skin
            %light; %add light for shading
%             
%             h = surf(obj.x_coord, obj.y_coord, obj.state);
%             view(az,el);
%             set(h,'LineStyle','none');
%             axis([0 obj.width 0 obj.length 0 obj.thickness*1.5]);
%           
            %C = cmap(end,:);            

%             fill3([0 0 obj.width obj.width],...
%                 [0 0 0 0],...
%                 [0 obj.thickness obj.thickness 0],...
%                 C);
%             
%             fill3([0 0 0 0],...
%                 [0 0 obj.length obj.length],...
%                 [0 obj.thickness obj.thickness 0],...
%                 C);
            %%%%% End Plot Tissue Sides
            
        end
        
        function laser = get_laser(obj)
            laser = obj.laser;
        end
        
        function r_th = calculate_max_ablation_r(obj)
            %returns r_th which is the radius of ablation effect of the 
            %laser for the largest Theta_0_th.
            %determine region of ablation surrounding laser
            %search for region where incident radiant exposure (Thet_0)
            %that is greater than Theta_th);
            
            r_th = 0; %threshold radius within which ablation occurs
            vox_area = (1/obj.resolution)^2; %mm^2
            %obj.Theta_th = obj.density * obj.h_abl / obj.mu_a;
            %Theta_th_min = min(min(obj.material_properties.density ...
            %                        .* obj.material_properties.h_abl ...
            %                        ./ obj.material_properties.mu_a));
            Theta_th_min = min(min(obj.irrad_thresh));
            Theta_0 = 2*Theta_th_min;
            r_temp = 0;
            while (Theta_0 > Theta_th_min)
                r_temp = r_temp + 1/obj.resolution;
                irrad_level = obj.laser.irradiance(r_temp,0);
                Theta_0 = irrad_level; %*dt/vox_area (I think)
            end
            
            r_th = r_temp;
        end
        
        %Static version now, can be updated to be more complex.
        function d = get_density(obj,x,y)
            d = obj.material_properties.density;
        end
        
        %Static version now, can be updated to be more complex.
        function h_abl = get_h_abl(obj,x,y)
            h_abl = obj.material_properties.h_abl;
        end
        
        %Static version now, can be updated to be more complex.
        function theta_th_xy = get_theta_th_xy(obj,x,y)
            theta_th_xy = obj.material_properties.theta_th_xy;
        end
        
        %minus: overloaded minus operator for Sample A - Sample B.
        %Subtracts State's. Must have same meshgrid (x_coord/y_coord)
        function aSample = minus(SampleA, SampleB)
            Zq = SampleA.state - SampleB.state;
            aSample = SampleA;
            aSample.state = Zq;
        end
    end %end methods section
end
