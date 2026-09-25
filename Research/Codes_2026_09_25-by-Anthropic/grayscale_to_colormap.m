function rgbImage = grayscale_to_colormap(grayImage, cmap)
    % Normalize grayscale image to [0, 1]
    grayImage = double(grayImage);
    grayImage = grayImage - min(grayImage(:));
    grayImage = grayImage / max(grayImage(:));

    % Number of colors in the colormap
    numColors = size(cmap, 1);

    % Scale grayscale image to indices in the colormap
    indices = round(grayImage * (numColors - 1)) + 1;

    % Initialize RGB image
    [rows, cols] = size(grayImage);
    rgbImage = zeros(rows, cols, 3);

    % Map grayscale to RGB using the custom colormap
    for i = 1:rows
        for j = 1:cols
            rgbImage(i,j,:) = cmap(indices(i,j), :);
        end
    end
end