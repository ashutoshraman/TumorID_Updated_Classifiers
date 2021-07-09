classdef laser405 < handle
    properties
        %width
        position
        %cutting_rate
        %beam_profile
        mode %continuous vs pulsed
        %pulse_length %
        power_level %0 to 10W for our CO2 Laser
        wattage %0 to 100%
        spot_size %Later this will be used.
        x_coord %indcices for beam profile
        y_coord %indicies for beam profile
        pwr_function %2D function describing intensity vs distance from beam center, usually Gaussian
        wavelength = 405; %not used for anything yet, just noting wavelength
    end
    
    methods
        function obj = laser405(wattage, spot_size)
            if exist('wattage','var')
                obj.wattage = wattage;
            else
                obj.wattage = 0.2; %watts
                disp('Laser: Setting duty cycle to 200ma (default).');
            end %duty cycle
            
            if exist('spot_size','var')
                obj.spot_size = spot_size;
            else
                obj.spot_size = 1.75;
                disp('Laser: Setting spot size to 1.75mm (default).');
            end %d
            
            %obj.spot_size = 0.5; %in mm
            %obj.power_level = 10; %0 to 10W for our CO2 Laser
            %obj.duty_cycle = 80; 
            obj.set_laser_settings(obj.wattage, obj.spot_size);
            %obj.width = obj.spot_size / 2;
            
            obj.position.x = 2; %this is initial center of beam in x (mm)
            obj.position.y = 2; %this is initial center of beam in y (mm)
            %obj.cutting_rate = 1; %dh/dt TODO: This is replaced when we include an ablation model (i.e. Beer Lambert)
            lspace = linspace(-3*obj.spot_size/2,3*obj.spot_size/2,99);
            [X,Y] = meshgrid(lspace,lspace);
            obj.x_coord = X;
            obj.y_coord = Y;
            
            %TODO: Make beam_profile a function of spot_size
            %obj.beam_profile = exp(-obj.power_level*((X.^2 + Y.^2));
            %obj.beam_profile = obj.power_level * exp(-1*((X/obj.width).^2+(Y/obj.width).^2));
        end
        
        function plot(obj)
            figure;
            lspace = linspace(-1.5*obj.spot_size,1.5*obj.spot_size,99);
            [X,Y] = meshgrid(lspace(1:2:end),lspace(1:2:end));
            %beam_profile = obj.power_level * exp(-1*((X/obj.spot_size).^2+(Y/obj.spot_size).^2));
            beam_profile = obj.pwr_function(sqrt((X.^2 + Y.^2)));
            surf(X, Y, beam_profile);
            shading interp;
            grid off;
            set(gcf, 'Color', 'White')
            xlabel('Width (mm)', 'FontSize', 15);
            ylabel('Length (mm)', 'FontSize', 15);
            zlabel('Fluence (W/mm^2)', 'FontSize', 15);
            title(['Peak Fluence = ' ...
                 num2str(obj.pwr_function(0))],'FontSize',15);
            
            %plot laser
            %[X Y Z] = cylinder(laser.width);
            %surf(X+laser.position.y, Y+laser.position.x,Z+sample.thickness*1.2);
            
        end
        
        function irrad = irradiance(obj, x, y)
            %irrad = obj.power_level * exp(-1*((x/obj.width).^2+(y/obj.width).^2));
            irrad = obj.pwr_function(sqrt(x.^2 + y.^2));
        end
        
        function set_power(obj, pwr_level)
            obj.power_level = pwr_level;
        end
        
        function set_laser_settings(obj, wattage, spot_size)
%             %%shortened code to quicken, remove later
%             obj.spot_size = spot_size;
%             obj.power_level = 9.1362;
%             w0 = spot_size/2;
%             obj.pwr_function = @(r) obj.power_level*exp(-(2*r.^2/w0^2)); 
%             return;
%             %%
            
            
            %disp('set_laser_settings');
            
            w0 = spot_size/2;
            obj.spot_size = spot_size;
            if (wattage <= 0 || wattage > 320)
                disp('ERROR: Wattage must be 0 < wattage < 320. Setting to 200.');
                wattage = 200;
            end
            %get peak power from duty cycle using function fit from data
                %collected from beam power experiment in Wes Hill's thesis
            %3rd degree polynomial fitted to mean of DC Power
            %pwr_level = @(x) -1.661e-6*(x^3)-0.0004575*x^2+0.1729*x +.9348;
            %total_power = pwr_level(wattage);
            
            total_power = wattage;
            
            %determine power distribution at that spot size
                %i.e. determine E0, peak intensity.
            %employ half-interval search for E0
            Q = 0.0;
            E0 = 0.0;
            low = 0;
            high = 320;
            epsilon = 0.01;
            while (abs(Q - total_power) > epsilon)
                E0 = (low+high)/2;
                pwr = @(r) E0*exp(-(2*r.^2/w0^2));
                pwr2 = @(x,y) pwr(sqrt(x.^2 + y.^2));
                %make sure that this equaiton above is correct and consistent with the other eqns.
                Q = integral2(pwr2, -inf, inf, -inf, inf);
                if Q > total_power
                    high = E0;
                elseif Q < total_power
                    low = E0;
                end
            end
            obj.power_level = E0;
            obj.pwr_function = @(r) obj.power_level*exp(-(2*r.^2/w0^2)); 
        end

    end %end methods section
end
