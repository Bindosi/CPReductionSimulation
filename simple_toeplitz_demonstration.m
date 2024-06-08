% Generate transmitted OFDM symbol (without CP)
num_carriers = 64;
num_taps = 7;
tx_symbols = randn(num_carriers, 1) + 1i * randn(num_carriers, 1);

% Generate channel impulse response
channel_response = randn(num_taps, 1) + 1i * randn(num_taps, 1);

% Add cyclic prefix
cp_length = 10; % Length of cyclic prefix
tx_symbols_cp = [tx_symbols(end-cp_length+1:end); tx_symbols];

% Construct Toeplitz matrix for channel convolution
H_toeplitz = toeplitz([channel_response; zeros(num_carriers+cp_length-1,1)], [channel_response(1); zeros(num_carriers-1,1)]);

% Convolution with the channel using Toeplitz matrix multiplication
rx_symbols = H_toeplitz * tx_symbols_cp;

% Remove cyclic prefix from received symbols
rx_symbols = rx_symbols(cp_length+1:end);

