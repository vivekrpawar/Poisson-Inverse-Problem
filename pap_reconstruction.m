image = imread('barbara256.png');
double(image);
% Get the size of the image
[height, width, ~] = size(image);   

beta = 0.5;
lambda = 0.25;
sigma = sqrt(beta/lambda); 

peak_value_set = [0.8];
noisy_image_snr = [];
recon_image_snr = [];

for i=1:length(peak_value_set)
    disp(i); 
    noisy_image = poisson_contaminated_image(image, peak_value_set(i));  
    reconstructed_image = PIP_reconstruction(imresize(noisy_image, 0.5), lambda, sigma);
    reconstructed_image = imresize(reconstructed_image, 2);
    snr = calculate_snr(image, noisy_image);
    recon_snr = calculate_snr(image, reconstructed_image);
    disp(["peak value", peak_value_set(i)]);
    disp(["snr", snr]);
    disp(["recon snr", recon_snr]);
    noisy_image_snr = [noisy_image_snr, snr];
    recon_image_snr = [recon_image_snr, recon_snr];
    rmse1 = norm(double(image(:))-double(noisy_image(:)))/norm(double(image(:)));
    rmse2 = norm(double(image(:))-double(reconstructed_image(:)))/norm(double(image(:)));
    disp(["rmse1", rmse1]);
    disp(["rmse2", rmse2]);
    figure;
    subplot(1,3,1);
    imshow(image, []);
    subplot(1,3,2);
    imshow(noisy_image, []);
    subplot(1,3,3);
    imshow(reconstructed_image, []);
end

%{
plot(peak_value_set, noisy_image_snr, 'b-', 'LineWidth', 2); % Blue solid line with a linewidth of 2
hold on; % Hold the plot to add another line

% Plot the second line
plot(peak_value_set, recon_image_snr, 'r--', 'LineWidth', 2); % Red dashed line with a linewidth of 2

% Add labels and legend
xlabel('Peak Value');
ylabel('SNR');
title('SNR value of noisy and reconstructed images (Cameraman)');
legend('Noisy Image SNR', 'Reconstructed Image SNR');

% Turn off the hold to allow further plot modifications
hold off;
%}

function reconstructed_image = PIP_reconstruction(y, lambda, sigma)
    [h, w] = size(y);
    y = double(y);
    x = double(zeros(h, w));
    u = double(zeros(h, w));
    v = double(ones(h, w))*(4*(sqrt(3/8)+1));
    i = double(ones(h, w));
    for i=1:50
        disp(i);
        x = ((lambda*(v-u-i) +  sqrt((lambda*(v-u)-i).^2 + 4*lambda*y))) / (2*lambda);  
        v = BM3D(x+u, sigma);
        u = u + (x - v);
    end
    reconstructed_image = x;
end
