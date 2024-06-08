clear
clc;
% Define system parameters
num_antennas_bs = 32;   % Number of base station antennas
num_carriers = 64;      % Number of OFDM carriers
num_taps = 7;           % Number of taps in the channel
num_symbols = 1;    % Number of OFDM symbols
SNR_dB = 0:5:30;        % Range of SNR values to simulate

% Cyclic prefix sizes to compare
cp_sizes = [15, 9];

% Initialize BER arrays for each CP size
BER = zeros(length(cp_sizes), length(SNR_dB));

% Simulation loop for different SNR values
for idx_cp = 1:length(cp_sizes)
    cp_size = cp_sizes(idx_cp);
    
    for idx_snr = 1:length(SNR_dB)
        % Generate random bits
        tx_bits = randi([0, 1], num_carriers * num_symbols, 1);
        
        % Modulation (BPSK)
        tx_symbols = 2 * tx_bits - 1;
        
        % Reshape symbols into symbol blocks
        tx_symbols_block = reshape(tx_symbols, num_carriers, num_symbols);
        
        % Generate channel matrix (num_antennas_bs x num_taps)
        H = (randn(num_antennas_bs, num_taps) + 1i * randn(num_antennas_bs, num_taps)) / sqrt(2);
        
        % Generate additive white Gaussian noise
        noise_power = 10^(-SNR_dB(idx_snr) / 10);
        noise = sqrt(noise_power) * (randn(num_antennas_bs, num_carriers + cp_size, num_symbols) + 1i * randn(num_antennas_bs, num_carriers + cp_size, num_symbols));
        
        % Generate cyclic prefix
        tx_symbols_cp = [tx_symbols_block(end-cp_size+1:end, :); tx_symbols_block];
        
        % Transmit through the channel
        rx_signal = sum(H * tx_symbols_cp(:,:,1), 1) + noise;
        
        % Remove cyclic prefix
        rx_signal_no_cp = rx_signal(:, cp_size+1:end, :);
        
        % Perform MRC combining in frequency domain
        W = conj(H) ./ norm(H, 'fro');
        combined_signal = sum(W .* rx_signal_no_cp, 1);
        
        % Demodulation (BPSK)
        rx_symbols = sign(real(combined_signal));
        
        % BER computation
        errors = sum(abs(tx_symbols_block - rx_symbols) > 0, 'all');
        BER(idx_cp, idx_snr) = errors / (num_carriers * num_symbols);
    end
end

% Plot BER vs SNR for each CP size
figure;
for idx_cp = 1:length(cp_sizes)
    semilogy(SNR_dB, BER(idx_cp, :), 'o-', 'LineWidth', 2);
    hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('CP Size 15', 'CP Size 9');
title('BER vs SNR for Different Cyclic Prefix Sizes');
