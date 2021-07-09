%This class holds a path to be used during a cut with the TumorCNC
%It contains a set of points to visit with their associated
%cut time and cut power. It also contains methods to manipulate
%the data, such as shifting it from the origin and rotating.
classdef cut_path_obj < matlab.mixin.Copyable
    properties
        location  = [0,0]; %(x,y) position of bottom left most point
        t_points %time at each point
        x_points %
        y_points %
        p_points %power
        currentPoint = 1;%between 1 and length(x_points), used to iterate through points
    end
    
    methods
        function obj = cut_path_obj(type, varargin)
            switch (nargin)
                case 2
                    default_args = varargin{1};
                case 1
                    disp('cut_path_obj: using default variables for path.');
                    L = 10; %length (mm)
                    W = 10; %width (mm)
                    d = 1.75; %beam width (mm)
                    alpha = 0.5; %step size in terms of r between raster cuts
                    dir = 'x'; %direction of cuts 'x' or 'y';
                    int_time = 5; %speed of cut mm/s (distance of cut path/speed gives int time) (do in terms of int time rather than speed)
                    %dt = 0.01; %Note: Maybe remove this as an arg??
                    step_size = 11; %points per centimeter
                    %default_args = {L,W,d,alpha,dir,speed,dt};
                    default_args = {L,W,d,alpha,dir,int_time,step_size};
            end
            
            cut_path = generate_path(type,default_args); 
            obj.t_points = cut_path(1,:);
            obj.x_points = cut_path(2,:);
            obj.y_points = cut_path(3,:);
            obj.p_points = 80*ones(1,size(cut_path,2)); %% why 80 here
        end
        
        function plot(obj)
           %This should visualize the cut points in a scatter plot
           %It would be nice to be able to show power and time/speed
           %of the cut some how. Maybe show power in color, and 
           %the speed of cut in terms of points-per-inch or something.
           figure;
           scatter(obj.x_points + obj.location(1), obj.y_points+obj.location(2), 12, obj.t_points, 'filled'); hold on;
           
           %plot start and end points
           scatter(obj.x_points(1) + obj.location(1), obj.y_points(1)+obj.location(2), 25, 'green','x');
           scatter(obj.x_points(end) + obj.location(1), obj.y_points(end)+obj.location(2), 25, 'red', 'x');

           %plot direction arrows
           dx = diff(obj.x_points);
           dy = diff(obj.y_points);
           quiver(obj.x_points(1:50:end) + obj.location(1), obj.y_points(1:50:end)+obj.location(2), dx(1:50:end), dy(1:50:end), ...
                'linewidth',1, 'color', 'blue');
           
           xlabel('Width [mm] (X)','FontSize',13); 
           ylabel('Length [mm] (Y)','FontSize',13);
           axis square;
        end
        
        function set_location(obj,xy)
           %xy = [x,y];
            x_median = (max(obj.x_points) + min(obj.x_points)) / 2;
            y_median = (max(obj.y_points) + min(obj.y_points)) / 2;
            obj.location = [xy - [x_median, y_median]];
        end
        
        function set_x_points(obj, x_pts)
            obj.x_points = reshape(x_pts, [1,length(x_pts)]);
            obj.check_values();
        end
            
        function set_y_points(obj, y_pts)
            obj.y_points = reshape(y_pts, [1,length(y_pts)]);
            obj.check_values();
        end
        
        function set_p_points(obj, p_pts)   
            obj.p_points = reshape(p_pts, [1,length(p_pts)]);
            obj.check_values();
        end
        
        function set_t_points(obj, t_pts)
            obj.t_points = reshape(t_pts, [1,length(t_pts)]);
            obj.check_values();
        end 
        
        function set_xy_points(obj, pts, speed) %change speed to path/int_time
            obj.x_points = pts(:,1);
            obj.y_points = pts(:,2);
    
            %caculate t points.
            distance_travled = cumsum(sqrt((obj.x_points(2:end)-obj.x_points(1:end-1)).^2 ...
                + (obj.y_points(2:end)-obj.y_points(1:end-1)).^2));
            T = distance_travled / speed;
            obj.t_points = [0; T];
            disp('Warning: Set the p_points too!');
        end
        
        function BC = export_backward_compatible(obj)
            %This function will be removed ASAP. three_d_laser_ablation()
            %needs to be setup to work with this new cut_path_obj
            %Note: Needs to be updated, the location doesn't work out.
            BC = [obj.t_points; ...
                obj.x_points + obj.location(1); ...
                obj.y_points + obj.location(2)];
        end
        
        function export_csv(obj, filename)
            center_at_origin = input('Center at origin? 0 No, 1 Yes');
            
            if center_at_origin == 0
                data = [reshape(obj.t_points, length(obj.t_points),1), ...
                    reshape(obj.x_points + obj.location(1), length(obj.x_points),1), ...
                    reshape(obj.y_points + obj.location(2), length(obj.y_points),1), ...
                    reshape(obj.p_points, length(obj.p_points),1)];
            else 
                data = [reshape(obj.t_points, length(obj.t_points),1), ...
                    reshape(obj.x_points, length(obj.x_points),1), ...
                    reshape(obj.y_points, length(obj.y_points),1), ...
                    reshape(obj.p_points, length(obj.p_points),1)];
            end
                
            csvwrite(filename, data);
        end
        
        function export_xyz_csv(obj, filename, z)
            center_at_origin = input('Center at origin? 0 No, 1 Yes');
            z_points = z * ones(size(obj.x_points));
            if center_at_origin == 0
                data = [...
                    reshape(obj.x_points + obj.location(1), length(obj.x_points),1), ...
                    reshape(obj.y_points + obj.location(2), length(obj.y_points),1), ...
                    reshape(z_points, length(obj.p_points),1)];
            else 
                data = [...
                    reshape(obj.x_points, length(obj.x_points),1), ...
                    reshape(obj.y_points, length(obj.y_points),1), ...
                    reshape(zeros(size(obj.x_points)), length(obj.x_points),1)];
            end
                
            csvwrite(filename, data);
         end
        
        %This function crops the cut points to be within a boundary
        function [in, on] = cropToBoundary(obj, boundary, speed)
            %Boundary = [x, y] (vectors of points)
            [in, on] = inpolygon(obj.x_points+obj.location(1), obj.y_points+obj.location(2),boundary(:,1),boundary(:,2));
            obj.x_points = obj.x_points(in | on);
            obj.y_points = obj.y_points(in | on);
            obj.p_points = obj.p_points(in | on);
            
            %%Note: this should be fixed / implemente
%             %caculate t points.
%             distance_travled = cumsum(sqrt((obj.x_points(2:end)-obj.x_points(1:end-1)).^2 ...
%                 + (obj.y_points(2:end)-obj.y_points(1:end-1)).^2));
%             T = distance_travled / speed;
%             obj.t_points = [0; T'];
            
            %Hack to make it work
            dt = obj.t_points(2) - obj.t_points(1);
            obj.t_points = cumsum(dt*ones(size(obj.p_points)));
        end
        
        function txy = nextPoint(obj)
            try
                txy = [obj.t_points(obj.currentPoint), ...
                    obj.x_points(obj.currentPoint) + obj.location(1), ...
                    obj.y_points(obj.currentPoint) + obj.location(2)];
                obj.currentPoint = obj.currentPoint + 1;
            catch
                disp('cut_path_obj: End of path reached. Resetting count. No point returned.');
                txy = NaN; %error, we're out of pts
                obj.currentPoint = 1;
            end
        end
        
        
        function txy = getPoint(obj,k)
            txy = [obj.t_points(k), ...
                    obj.x_points(k) + obj.location(1), ...
                    obj.y_points(k) + obj.location(2)];
        end
        
        %this function checks the values to make sure they're correct
        function check_values(obj)
            L = length(obj.x_points);
            if (        (length(obj.y_points) ~= L) ...
                    || (length(obj.t_points) ~= L) ...
                    || (length(obj.p_points) ~= L) )
                disp('Error! cut_path_obj.check_values() failed. Check t,x,y,p.');
            end
        end
        
    end %end methods section
end