function snr_value = calculate_snr(clean_image, noisy_image)
    % Convert images to double for calculations
    clean_image = double(clean_image);
    noisy_image = double(noisy_image);
    
    % Calculate signal power (using clean image)
    signal_power = sum(clean_image(:) .^ 2) / numel(clean_image);
    
    % Calculate noise power (using the difference between noisy image and clean image)
    noise = noisy_image - clean_image;
    noise_power = sum(noise(:) .^ 2) / numel(clean_image);
    
    % Calculate SNR
    snr_value = 10 * log10(signal_power / noise_power);
end
