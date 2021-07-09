% % Create a logical image of a circle with specified
% % diameter, center, and image size.
% % First create the image.
% imageSizeX = 640;
% imageSizeY = 480;
% [columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% % Next create the circle in the image.
% centerX = 320;
% centerY = 240;
% radius = 100;
% circlePixels = (rowsInImage - centerY).^2 ...
%     + (columnsInImage - centerX).^2 <= radius.^2;
% % circlePixels is a 2D "logical" array.
% % Now, display it.
% image(circlePixels) ;
% colormap([0 0 0; 1 1 1]);
% title('Binary image of a circle');

function custom_lesion = create_lesion(old_state, length, width, thickness, shape)
    [x_dim y_dim] = size(old_state);
    [columnsInImage rowsInImage] = meshgrid(1:y_dim, 1:x_dim);
    if shape == 'circle'
        centerX = width/2;
        centerY = length/2;
        if length == width
            radius = length/4;
            mask = (rowsInImage - centerX).^2 ...
                + (columnsInImage - centerY).^2 <= radius.^2;
        else
            radiusX = width/4;
            radiusY = length/4;
            mask = (rowsInImage - centerX).^2 ./ radiusX^2 ...
                + (columnsInImage - centerY).^2 ./ radiusY^2 <= 1;
        end
        custom_lesion = old_state .* mask;
        custom_lesion(custom_lesion == thickness) = thickness/2;
        custom_lesion(custom_lesion == 0) = thickness;
    elseif shape == 'square'
        old_state(length/4:length/4 + length/2, width/4:width/4 + width/2) = thickness / 2;
        custom_lesion = old_state;
    end
end