function [blurred_image, H, blurKernel] = image_with_uniform_noise(image)
    [height, width] = size(image);
    disp([height, width]);
    % Define the size of the image and kernel
    % Define the patch size
    patch_size = 32;

    kernelSize = 9; % Adjust according to the size of your blur kernel
    
    % Create a blur kernel (you can define your own blur kernel)
    blurKernel = fspecial('average', [kernelSize, kernelSize]); % Uniform blur kernel 
    %blurKernel = circshift(blurKernel, -3, 2);
    blur_vec_matrix = zeros(32);

    blur_vec_matrix(1:kernelSize, 1:kernelSize) = blurKernel; 
    blur_matrix_row = blur_vec_matrix(:); 
    % blur_matrix_row = circshift(blur_matrix_row, -300, 1);
    % Fill the circulant matrix
    disp(sum(blur_matrix_row));
    n = patch_size*patch_size; 
    circMatrix = zeros(n);
    for i = 1:n
        circMatrix(:, i) = circshift(blur_matrix_row, i-1);
    end   
    
    H = circMatrix; 
    
    
    % Compute the number of patches in each dimension
    num_patches_vertical = floor(height / patch_size);
    num_patches_horizontal = floor(width / patch_size);   
    
    blurred_image = zeros(size(image));
    for i = 1:num_patches_vertical
        for j = 1:num_patches_horizontal
            % Compute the coordinates of the current patch
            row_start = (i - 1) * patch_size + 1;
            row_end = row_start + patch_size - 1;
            col_start = (j - 1) * patch_size + 1;
            col_end = col_start + patch_size - 1;
            
            % Extract the patch
            patch = image(row_start:row_end, col_start:col_end, :);
    
            % Store the patch in the cell array
            vectorized_patch = patch(:);  
            blurred_vector = H * double(vectorized_patch);
    
            % Insert the patch into the reconstructed image
            blurred_patch = reshape(blurred_vector, 32, 32);
            blurred_image(row_start:row_end, col_start:col_end, :) = blurred_patch;
        end
    end
end
