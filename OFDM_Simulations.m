clear; clc;
deltaF        = 15000;
Mod_Order     = 2;           % BPSK (1), QPSK(2), 8-PSK(3), 16-QAM(4), 64-QAM(6)         
N_fft         = 64;          % FFT size
Symbol_Period = 1/(deltaF*N_fft);

NumSym    = 10;        % no of realization, how many ofdm symbols do we transmit

SNR                 = (0:2:30);
BER_n_1             = zeros(1,length(SNR));
BER_n_2             = zeros(1,length(SNR));
BER_p_1             = zeros(1,length(SNR));
BER_p_2             = zeros(1,length(SNR));
BER_nc              = zeros(1,length(SNR));
BER_pc              = zeros(1,length(SNR));
Optimize_Proposed   = zeros(1,length(SNR));
Optimize_Normal     = zeros(1,length(SNR));

% Channel
% Multitap channel
PowerdBCluser1  =  [-1 -7 -17 -21 -25];  % Cluster 1 tap power profile 'dB'
DelayCluster1   =  [0 2 5 6 8];          % Cluster 2 delay 'sample'
PowerdBCluser2  =  [-3 -9 -19 -25 -31];  % Cluster 1 tap power profile 'dB'
DelayCluster2   =  [8 9 11 13 14];       % Cluster 2 delay 'sample'
VR1_Antennas    =   32;
VR2_Antennas    =   22;

[H_cir_1, Lch_1, H_cir_2, Lch_2, Lch_Combined]    =  Channel(PowerdBCluser1,PowerdBCluser2,DelayCluster1,DelayCluster2,NumSym,VR1_Antennas,VR2_Antennas);

H_cfr_1           =  fft([H_cir_1; zeros(N_fft-Lch_1,NumSym,VR1_Antennas)]);
H_cfr_2           =  fft([H_cir_2; zeros(N_fft-Lch_2,NumSym,VR2_Antennas)]);

% Determining the CP's
CP_Normal    = Lch_Combined;
CP_Proposed  = max(Lch_1,Lch_2);
Lch_selected = CP_Proposed;

% Delay in time micro seconds
CP_time_normal      = CP_Normal     * 10^-8;
CP_time_proposed    = CP_Proposed   * 10^-8;

% Symbol size (duration)
N_sym              = N_fft + CP_Normal;   % Symbol duration
N_sym_proposed     = N_fft + CP_Proposed; % Symbol duration

% FFT Process
I           = eye(N_fft);  
F           = fft(I);     % FFT matrix
IF          = ifft(I);    % IFFT matrix
CP          = [I(N_fft-CP_Normal+1:N_fft,:) ; I];        % CP Insetion Matrix
CP_P        = [I(N_fft-CP_Proposed+1:N_fft,:) ; I];      % CP Insetion Matrix Proposed
T           = [ zeros(N_fft,CP_Normal) I];               % CP removal Matrix
T_proposed  = [ zeros(N_fft,CP_Proposed) I];             % CP removal Matrix Proposed

% Tx: Signal generation
X_bits    = randn(1,NumSym*N_fft*Mod_Order)>0;
X_modSyms = Modulator(X_bits, Mod_Order);
X_modSyms = reshape(X_modSyms,N_fft,NumSym);

% OFDM modulation
x             = sqrt(N_fft)*IF*X_modSyms;
x_cp_Normal   = CP*x;
x_cp_Proposed = CP_P*x;

% Pre allocation of matrices
y_proposed_1  =   zeros(N_sym_proposed, NumSym,VR1_Antennas);
y_proposed_2  =   zeros(N_sym_proposed, NumSym,VR2_Antennas);

y_normal_1    =   zeros(N_sym, NumSym,VR1_Antennas);
y_normal_2    =   zeros(N_sym, NumSym,VR2_Antennas);

% building the received symbols matrix for normal CP
   H_temp_1       =   zeros(N_sym + Lch_1 , N_sym, VR1_Antennas);  
   H_temp_1_p     =   zeros(N_sym_proposed + Lch_1, N_sym_proposed,VR1_Antennas);

    for sym       =   1:NumSym
        for ant   =   1:VR1_Antennas
            for s =   1:N_sym
                H_temp_1(s:s + Lch_1-1, s,ant)        = H_cir_1(:,sym,ant);
            end
             for sp =   1:N_sym_proposed
                H_temp_1_p(sp:sp + Lch_1-1,sp,ant)    = H_cir_1(:,sym,ant);
            end
        end
        H_toep_1                        = H_temp_1(1:N_sym,:,:);           % Channel Toeplitz Matrix
        H_toep_1_p                      = H_temp_1_p(1:N_sym_proposed,:,:);% Channel Toeplitz Matrix   
        
         
        y_normal_1(:,sym,ant)           = H_toep_1(:,:,ant)   * x_cp_Normal(:,sym); % Multiplying xcp sym's row and htoep creating TX OFDM symbols 
        y_proposed_1(:,sym,ant)         = H_toep_1_p(:,:,ant) * x_cp_Proposed(:,sym);
    end
   

% building the received symbols matrix for normal and proposed CP
H_temp_2        =   zeros(N_sym + Lch_2, N_sym,VR2_Antennas);
H_temp_2_p      =   zeros(N_sym_proposed + Lch_2, N_sym_proposed,VR2_Antennas);

for sym          =   1:NumSym
    for ant = 1:VR2_Antennas
            for samp =   1:N_sym
                H_temp_2(samp:samp + Lch_2-1,samp,ant)            = H_cir_2(:,sym,ant);
            end
             for sampp =   1:N_sym_proposed
                H_temp_2_p(sampp:sampp + Lch_2-1,sampp,ant)       = H_cir_2(:,sym,ant);
            end
    
        H_toep_2_p                   = H_temp_2_p(1:N_sym_proposed,:,:);               % Channel Toeplitz Matrix
        H_toep_2                     = H_temp_2(1:N_sym,:,:);                          % Channel Toeplitz Matrix
        y_proposed_2(:,sym,ant)      = H_toep_2_p(:,:,ant) * x_cp_Proposed(:,sym); 
        y_normal_2(:,sym,ant)        = H_toep_2(:,:,ant)   * x_cp_Normal(:,sym); % Multiplying xcp sym's row and htoep creating TX OFDM symbols
    end
end


Spectral_Efficiency_Normal      = Symbol_Period/(Symbol_Period  + CP_time_normal);
Spectral_Efficiency_Proposed    = Symbol_Period/(Symbol_Period  + CP_time_proposed );

mrcY_n_1     = zeros(N_fft, 1, VR1_Antennas);
mrcY_p_1     = zeros(N_fft, 1, VR1_Antennas);
mrcY_n_2     = zeros(N_fft, 1, VR2_Antennas);
mrcY_p_2     = zeros(N_fft, 1, VR2_Antennas);

for k=1:length(SNR)
    %AWGN
    y_p_rx_cp_1 = awgn(y_proposed_1,SNR(k),'measured');% Adding AWGN to the received OFDM symbols
    y_p_rx_cp_2 = awgn(y_proposed_2,SNR(k),'measured');% Adding AWGN to the received OFDM symbols

    y_n_rx_cp_1 = awgn(y_normal_1,SNR(k),'measured');  % Adding AWGN to the received OFDM symbols
    y_n_rx_cp_2 = awgn(y_normal_2,SNR(k),'measured');  % Adding AWGN to the received OFDM symbols
    
    %RX
    y_p_rx_1    = pagemtimes(repmat(T_proposed,1,1,VR1_Antennas),  y_p_rx_cp_1);     % CP removal
    y_p_rx_2    = pagemtimes(repmat(T_proposed,1,1,VR2_Antennas) , y_p_rx_cp_2);    % CP removal
    y_n_rx_1    = pagemtimes(repmat(T,1,1,VR1_Antennas) , y_n_rx_cp_1);             % CP removal
    y_n_rx_2    = pagemtimes(repmat(T,1,1,VR2_Antennas) , y_n_rx_cp_2);             % CP removal

    Y_p_1    = pagemtimes(repmat(F,1,1,VR1_Antennas),  y_p_rx_1/sqrt(N_fft)); % FFT processing
    Y_p_2    = pagemtimes(repmat(F,1,1,VR2_Antennas) , y_p_rx_2/sqrt(N_fft)); % FFT processing
    Y_n_1    = pagemtimes(repmat(F,1,1,VR1_Antennas) , y_n_rx_1/sqrt(N_fft)); % FFT processing
    Y_n_2    = pagemtimes(repmat(F,1,1,VR2_Antennas) , y_n_rx_2/sqrt(N_fft)); % FFT processing

    Norm_H_1            = abs(H_cfr_1).^2;
    Norm_H_1_sum        = sum(Norm_H_1,3);
    Norm_H_1_sum_3D     = repmat(Norm_H_1_sum,1,1,VR1_Antennas);

    Norm_H_2            = abs(H_cfr_2).^2;
    Norm_H_2_sum        = sum(Norm_H_2,3);
    Norm_H_2_sum_3D     = repmat(Norm_H_2_sum,1,1,VR2_Antennas);
    
    mrcY_n_1_3D         = Y_n_1.*conj(H_cfr_1)./Norm_H_1_sum_3D;
    mrcY_p_1_3D         = Y_p_1.*conj(H_cfr_1)./Norm_H_1_sum_3D;
    mrcY_n_2_3D         = Y_n_2.*conj(H_cfr_2)./Norm_H_2_sum_3D;
    mrcY_p_2_3D         = Y_p_2.*conj(H_cfr_2)./Norm_H_2_sum_3D;

    mrcY_n_1            = sum(mrcY_n_1_3D,3);
    mrcY_p_1            = sum(mrcY_p_1_3D,3);
    mrcY_n_2            = sum(mrcY_n_2_3D,3);
    mrcY_p_2            = sum(mrcY_p_2_3D,3);


    Y_normal_combined   = mrcY_n_1 + mrcY_n_2;
    Y_proposed_combined = mrcY_p_1 + mrcY_p_2;

    Y_demodSyms_n_1     = reshape(mrcY_n_1,1,N_fft*NumSym);
    Y_demodSyms_p_1     = reshape(mrcY_p_1,1,N_fft*NumSym);
    Y_bits_n_1          = DeModulator(Y_demodSyms_n_1,Mod_Order);
    Y_bits_p_1          = DeModulator(Y_demodSyms_p_1,Mod_Order); 

    Y_demodSyms_nc      = reshape(Y_normal_combined,1,N_fft*NumSym);
    Y_demodSyms_pc      = reshape(Y_proposed_combined,1,N_fft*NumSym);

    Y_demodSyms_n_2     = reshape(mrcY_n_2,1,N_fft*NumSym);
    Y_demodSyms_p_2     = reshape(mrcY_p_2,1,N_fft*NumSym);
    Y_bits_n_2          = DeModulator(Y_demodSyms_n_2,Mod_Order); 
    Y_bits_p_2          = DeModulator(Y_demodSyms_p_2,Mod_Order); 
    
    Y_bits_nc           = DeModulator(Y_demodSyms_nc,Mod_Order); 
    Y_bits_pc           = DeModulator(Y_demodSyms_pc,Mod_Order); 
    % Optimization function
    Optimize_Normal(:,k)     = Spectral_Efficiency_Normal    * log (1+SNR(k));
    Optimize_Proposed(:,k)   = Spectral_Efficiency_Proposed  * log (1+SNR(k));
    %BER calculation
    [~,BER_n_1(k)]           = biterr(X_bits,Y_bits_n_1);
    [~,BER_n_2(k)]           = biterr(X_bits,Y_bits_n_2);
    [~,BER_p_1(k)]           = biterr(X_bits,Y_bits_p_1);
    [~,BER_p_2(k)]           = biterr(X_bits,Y_bits_p_2);

    [~,BER_nc(k)]            = biterr(X_bits,Y_bits_nc);
    [~,BER_pc(k)]            = biterr(X_bits,Y_bits_pc);
  
end

figure
    semilogy(SNR,BER_n_1,'-','Color','g','LineWidth',1.5)
    hold on
    semilogy(SNR,BER_n_2,'*','Color','g','LineWidth',1.5)
    hold on
    semilogy(SNR,BER_p_1,'-','Color','r','LineWidth',1.5)
    hold on
    semilogy(SNR,BER_p_2,'*','Color','r','LineWidth',1.5)
    hold on
    semilogy(SNR,BER_nc,'o','LineWidth',1.5)
    hold on
    semilogy(SNR,BER_pc,'--','LineWidth',1.5)
   
    legend('Normal VR1','Normal VR2','Proposed VR1','Proposed VR2', 'Combined Norm', 'Combined Prop','Location', 'SouthEast')
    xlabel('SNR (dB)'); ylabel('BER')
    grid on
    box on

figure
    semilogy(SNR,Optimize_Normal,'-','Color','g','LineWidth',1.5)
    hold on
    semilogy(SNR,Optimize_Proposed,'*','Color','r','LineWidth',1.5)
    
    legend('Normal','Proposed','Location', 'SouthEast')
    
    xlabel('SNR (dB)'); ylabel('Max function (Spec Eff * log(1+SNR)')
    grid on
    box on
