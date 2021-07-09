function cut_path = generate_path(type, arg)
    if strcmp(arg,'default')
        %set defaul args
        if strcmp(type,'raster')
            L = 5; %length (mm)
            W = 5; %width (mm)
            d = 1.75; %beam width (mm)
            alpha = 0.25; %step size in terms of r between raster cuts
            dir = 'y'; %direction of cuts 'x' or 'y';
            int_time = 5; %path length over integration time gives speed of cut
            step_size = 50; %points per cm
            arg = {L,W,d,alpha,dir,int_time,step_size};
        elseif strcmp(type,'circle')
            r = 3; %Radius in mm
            step_size = 50; %points per cm
            speed = 5; %speed of cut mm/s %change to int_time
            arg = {r,step_size, speed};
        end
    end
    
    if strcmp(type, 'circle')
        cut_path = generate_circle_path(arg); %arg = {r, step_size, speed}
    elseif strcmp(type, 'raster')
        cut_path = generate_raster_path(arg);
    elseif strcmp(type,'points')
        cut_path = [arg{1},arg{2}, arg{3}]'; %arg = {t,X,Y}
    elseif strcmp(type,'empty')
        cut_path = generate_empty_path();
    else
        %default to circlce for now.
        r = arg;
        cut_path = generate_circle_path(r);
    end
end

function cut_path = generate_empty_path(arg)
    disp('generate_path.generate empty path. Make this!');
    disp('Generating empty cut path.');
    cut_path = [nan; nan; nan]; %[t,X,Y]
end

function cut_path = generate_raster_path(arg)
    L = arg{1}; %length (mm)
    W = arg{2}; %width (mm)
    d = arg{3}; %beam width (mm)
    alpha = arg{4}; %step size in terms of r between raster cuts
    dir = arg{5}; %direction of cuts 'x' or 'y';
    int_time = arg{6}; 
    step_size = arg{7};
    if dir == 'x'
        num_cuts = floor(L / (alpha * d))+1; 
   
        if num_cuts == 1
            y_coords = 0;
        else
            y_coords = linspace(-L/2, L/2, num_cuts);
        end
        
        x_coords = linspace(-W/2, W/2, W * (1 / step_size));
        
        [X,Y] = meshgrid(x_coords, y_coords);
        X_temp = X';
        Y_temp = Y';
        
        %X = reshape(X',1,numel(X));
        %Y = reshape(Y',1,numel(Y));
        %dt = (X(2) - X(1)) / speed;
        
        %add code here to make a zig-zag
        X1 = X_temp; X1(:,1:2:end) = 0;
        X2 = X_temp; X2(:,2:2:end) = 0;
        Xnew = X1 + flip(X2);

        Y1 = Y_temp; Y1(:,1:2:end) = 0;
        Y2 = Y_temp; Y2(:,2:2:end) = 0;
        Ynew = Y1 + flip(Y2);   
        
        speed = W / int_time;
        
        
    else %dir == 'y'
        num_cuts = floor(W / (alpha * d))+1;
        y_coords = linspace(-L/2, L/2, L * (1 / step_size));
        
        if num_cuts == 1
            x_coords = 0;
        else
            x_coords = linspace(-W/2, W/2, num_cuts);
        end
        
        [X,Y] = meshgrid(x_coords, y_coords);
        X_temp = X;
        Y_temp = Y;
      
        %X = reshape(X,1,numel(X))
        %Y = reshape(Y,1,numel(Y));
        %dt = (Y(2) - Y(1)) / speed;
        
        %add code here to make a zig-zag
        X1 = X_temp; X1(:,1:2:end) = 0;
        X2 = X_temp; X2(:,2:2:end) = 0;
        Xnew = X1 + flip(X2);

        Y1 = Y_temp; Y1(:,1:2:end) = 0;
        Y2 = Y_temp; Y2(:,2:2:end) = 0;
        Ynew = Y1 + flip(Y2);
        
        speed = L / int_time;
    end
   
    
    Xnew = flip(reshape(Xnew,1,numel(Xnew)));
    Ynew = flip(reshape(Ynew,1,numel(Ynew)));
        
%     %%Old Code
%     %dt = 1 / (speed * (1/step_size) / 10);   %speed is (mm/s) pp_cc is points/cm
%     t = cumsum(dt * ones(1,length(X)));
%     cut_path = [t; X; Y];

    %%New Code
    dt = norm([Xnew(2)-Xnew(1); Ynew(2) - Ynew(1)]) / speed;
    t = cumsum(dt * ones(1,length(Xnew)));
    cut_path = [t; Xnew; Ynew];
end

function cut_path = generate_circle_path(arg) %see comments
    %returns an 2 by N array of (x;y) positions
    %arg = {r, step_size, speed}
    r = arg{1};
    step_size = arg{2};
    speed = arg{3}; % change to int_time
    
    distance = 2 * pi * r; %mm
    num_points = distance * (1 / step_size);
    arc_length = distance / num_points;
    ang = linspace(0,2*pi,2*pi*r/arc_length);
    xp = r*cos(ang(1:end-1));
    yp = r*sin(ang(1:end-1));
    
    dt = 1 / (speed * (1 / step_size) / 10);   %speed is (mm/s) pp_cc is points/cm  % change to speed=(2*pi*r/int_time)?
    t = cumsum(dt * ones(1,length(xp)));
    cut_path = [t; xp; yp];
end