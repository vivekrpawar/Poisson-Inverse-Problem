image = imread("barbara256.png");
disp(size(image));
%image = rgb2gray(image);
[height, width, ~] = size(image);
peak_value_set = [1];
original_image_snr = [];
reconstructed_image_snr = [];
for k = 1:length(peak_value_set)  
    [blurred_image, H, ~] = image_with_gaussian_noise(image);
    [snr1, snr2] = reconstruct_image(image, peak_value_set(k), height, width, blurred_image, H);
    disp(["peak value:", peak_value_set(k)]);
    disp(["Noisy image snr:", snr1]);
    disp(["Reconstructed image snr: ", snr2]);
    original_image_snr = [original_image_snr, snr1];
    reconstructed_image_snr = [reconstructed_image_snr, snr2];
end


plot(peak_value_set, original_image_snr, 'b-', 'LineWidth', 2); % Blue solid line with a linewidth of 2
hold on; % Hold the plot to add another line

% Plot the second line
plot(peak_value_set, reconstructed_image_snr, 'r--', 'LineWidth', 2); % Red dashed line with a linewidth of 2

% Add labels and legend
xlabel('Peak Value');
ylabel('SNR');
title('SNR value of noisy and reconstructed images (Cameramen) without binning.');
legend('Noisy Image SNR', 'Reconstructed Image SNR');

% Turn off the hold to allow further plot modifications
hold off;

function [snr1, snr2] = reconstruct_image(image, peak_value, height, width, blurred_image, H) 
    noisy_image = poisson_contaminated_image(blurred_image, peak_value); 
    patch_size = 32;
    beta = 1;
    lambda = 1;
    sigma = sqrt(beta/lambda);
    %h = imresize(H, 0.25);
    
    % Compute the number of patches in each dimension
    num_patches_vertical = floor(height / patch_size);
    num_patches_horizontal = floor(width / patch_size);
    reconstructed_image = zeros(size(image));
    num_of_iterations = 5;
    for i = 1:num_patches_vertical
        disp(i);
        for j = 1:num_patches_horizontal
            % Compute the coordinates of the current patch
            row_start = (i - 1) * patch_size + 1;
            row_end = row_start + patch_size - 1;
            col_start = (j - 1) * patch_size + 1;
            col_end = col_start + patch_size - 1;
            
            % Extract the patch
            patch = noisy_image(row_start:row_end, col_start:col_end, :); 
            %patch = imresize(patch, 0.5);
            % Store the patch in the cell array
   
            reconstructed_patch = zeros(size(patch));
            for k = 1: num_of_iterations
                reconstructed_patch = reconstructed_patch + PIP_reconstruction(patch, lambda, sigma, H);
            end
            reconstructed_patch = reconstructed_patch / num_of_iterations;
            %reconstructed_patch = imresize(reconstructed_patch, 2);
            % Insert the patch into the reconstructed image 
            reconstructed_image(row_start:row_end, col_start:col_end, :) = reconstructed_patch;
        end
    end
    snr1 = calculate_snr(image, noisy_image);
    snr2 = calculate_snr(image, reconstructed_image);
end

function reconstructed_image = PIP_reconstruction(y, lambda, sigma, H)
    lr=0.1;
    [h, w] = size(y);
    y = double(y);
    x = y;
    u = double(zeros(h, w));
    v = x;
    I = double(ones(h, w));
    I = I(:);
    for i=1:10
        %{
        for j=1:10 
            L= (-1) * H'*(y./(H*x)) + H'*I + lambda*(x-v+u); 
            x = x + lr*L; 
        end
        %}
        y = y(:);
        x = x(:);
        u = u(:);
        v = v(:);  
        lb = -inf*ones(size(x));
        ub = inf*ones(size(x));
        [n, ~] = size(x);
        % Define initial guess 
        
        % Call L-BFGS solver
        [x,f,info] = lbfgsb(@(x) objective_function(x, y, H, v, u, lambda), lb, ub); 
        % objective = @(x) (- y' * log(H * x) + ones(1, length(x)) * H * x + lambda / 2 * (x - v + u)' * (x - v + u));
        % x0 = x;
        % Solve the optimization problem using fminunc
        % options = optimoptions('fminunc', 'Display', 'iter'); % Adjust display option as needed
        % [x_opt, fval] = fminunc(objective, x0, options);
        y = reshape(y, h, w);
        x = reshape(x, h, w);
        u = reshape(u, h, w); 
        v = BM3D(x+u, sigma);
        u = u + (x - v);
    end
    reconstructed_image = x;
end


function [f, grad] = objective_function(x, y, H, v_k, u_k, lambda)
    % Compute the objective value
    f = - y' * log(H * x) + ones(1, length(x)) * H * x + (lambda / 2) * (x - v_k + u_k)' * (x - v_k + u_k);
    
    % Compute the gradient
    grad = - H' * (y ./ (H * x)) + H' * ones(length(x), 1) + lambda * (x - v_k + u_k); 
end