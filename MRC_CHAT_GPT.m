% Define system parameters
num_antennas_bs = 32;   % Number of base station antennas
num_carriers = 64;      % Number of OFDM carriers
num_taps = 7;           % Number of taps in the channel
num_symbols = 10000;    % Number of OFDM symbols
SNR_dB = 0:5:30;        % Range of SNR values to simulate

% Initialize BER arrays
BER_cp15 = zeros(size(SNR_dB));
BER_cp9 = zeros(size(SNR_dB));

% Simulation loop for different SNR values
for idx_snr = 1:length(SNR_dB)
    % Generate random bits
    tx_bits = randi([0, 1], num_carriers * num_symbols, 1);
    
    % Modulation (BPSK)
    tx_symbols = 2 * tx_bits - 1;
    
    % Reshape symbols into symbol blocks
    tx_symbols_block = reshape(tx_symbols, num_carriers, num_symbols);
    
    % Generate channel matrix (num_antennas_bs x num_taps)
    H = (randn(num_antennas_bs, num_taps) + 1i * randn(num_antennas_bs, num_taps)) / sqrt(2);
    
    % Generate cyclic prefix (CP) of size 15
    cp_size_15 = 15;
    tx_signal_cp15 = zeros(num_antennas_bs, num_carriers + cp_size_15, num_symbols);
    for idx_symbol = 1:num_symbols
        tx_signal_cp15(:, :, idx_symbol) = [tx_symbols_block(:, idx_symbol(end - cp_size_15 + 1:end)), tx_symbols_block(:, :, idx_symbol)];
    end
    
    % Generate cyclic prefix (CP) of size 9
    cp_size_9 = 9;
    tx_signal_cp9 = zeros(num_antennas_bs, num_carriers + cp_size_9, num_symbols);
    for idx_symbol = 1:num_symbols
        tx_signal_cp9(:, :, idx_symbol) = [tx_symbols_block(:, idx_symbol(end - cp_size_9 + 1:end)), tx_symbols_block(:, :, idx_symbol)];
    end
    
    % Transmit through the channel
    rx_signal_cp15 = zeros(num_antennas_bs, num_carriers + cp_size_15, num_symbols);
    rx_signal_cp9 = zeros(num_antennas_bs, num_carriers + cp_size_9, num_symbols);
    for idx_symbol = 1:num_symbols
        % Additive White Gaussian Noise (AWGN) with specified SNR
        noise = 1/sqrt(2) * (randn(num_antennas_bs, num_carriers + cp_size_15) + 1i * randn(num_antennas_bs, num_carriers + cp_size_15));
        noise_power = 10^(-SNR_dB(idx_snr) / 10);
        noise = sqrt(noise_power) * noise;
        
        % Received signal with CP of size 15
        rx_signal_cp15(:, :, idx_symbol) = sum(H * squeeze(tx_signal_cp15(:, :, idx_symbol)), 2) + noise;
        
        % Received signal with CP of size 9
        rx_signal_cp9(:, :, idx_symbol) = sum(H * squeeze(tx_signal_cp9(:, :, idx_symbol)), 2) + noise;
    end
    
    % Remove cyclic prefix
    rx_signal_cp15_no_cp = rx_signal_cp15(:, cp_size_15 + 1:end, :);
    rx_signal_cp9_no_cp = rx_signal_cp9(:, cp_size_9 + 1:end, :);
    
    % Perform MRC combining
    W = H' / norm(H, 'fro');
    combined_signal_cp15 = sum(W * rx_signal_cp15_no_cp, 1);
    combined_signal_cp9 = sum(W * rx_signal_cp9_no_cp, 1);
    
    % Demodulation (BPSK)
    rx_symbols_cp15 = sign(real(combined_signal_cp15));
    rx_symbols_cp9 = sign(real(combined_signal_cp9));
    
    % BER computation
    errors_cp15 = sum(abs(tx_symbols_block(:) - rx_symbols_cp15(:)) > 0);
    errors_cp9 = sum(abs(tx_symbols_block(:) - rx_symbols_cp9(:)) > 0);
    BER_cp15(idx_snr) = errors_cp15 / (num_carriers * num_symbols);
    BER_cp9(idx_snr) = errors_cp9 / (num_carriers * num_symbols);
end

% Plot BER vs SNR
figure;
semilogy(SNR_dB, BER_cp15, 'bo-', 'LineWidth', 2);
hold on;
semilogy(SNR_dB, BER_cp9, 'rs-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('CP Size 15', 'CP Size 9');
title('BER vs SNR for Different Cyclic Prefix Sizes');
