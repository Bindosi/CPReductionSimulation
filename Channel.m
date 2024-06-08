function [H_cir_1, Lch_1, H_cir_2, Lch_2, Lch_Combined] = Channel(PowerdBCluser1,PowerdBCluser2,DelayCluster1,DelayCluster2,N_realization,VR1_Antennas,VR2_Antennas)
% clear; clc;
% PowerdBCluser1  =  [-2 -8 -17 -21 -25]; % Channel tap power profile 'dB'
% DelayCluster1    =  [0 3 5 6 8];         % Channel delay 'sample'
% PowerdBCluser2  =  [-3 -9 -19 -25 -31];  % Cluster 1 tap power profile 'dB'
% DelayCluster2   =  [8 9 11 13 14];       % Cluster 2 delay 'sample'
% N_realization = 1000;

% ExcessDelay1    =  DelayCluster1(end) - DelayCluster1(1);
% ExcessDelay2    =  DelayCluster2(end) - DelayCluster2(1);
% SelectedExcessDelay   = max(ExcessDelay2,ExcessDelay1);
% x_Normal_CP     =   NormalTaps+1;
% x_Proposed_CP   =   ProposedTaps+1;


CombinedDelays    = [DelayCluster1 DelayCluster2];   % Combined delays of the two clusters
Time_adjusted_d2  = DelayCluster2 - DelayCluster2(1);

Power_1           =  10.^(PowerdBCluser1/10);
Power_2           =  10.^(PowerdBCluser2/10);

Ntap_1            =  length(PowerdBCluser1);
Ntap_2            =  length(PowerdBCluser2);

Lch_1             =  DelayCluster1(end) + 1;
Lch_2             =  Time_adjusted_d2(end) + 1;
Lch_Combined      =  CombinedDelays(end)+ 1;   % Channel length

h_1               =  zeros(Lch_1,1);
h_2               =  zeros(Lch_2,1);

H_cir_1           =  zeros(Lch_1,N_realization);
H_cir_2           =  zeros(Lch_2,N_realization); % either make this 15 or



 for ant = 1:VR1_Antennas
    %     rng(0,'twister')
    for Sym = 1:N_realization
        % creating channel in time domain complex entries using number of taps
        channel_1    =  (randn(1,Ntap_1)+ 1i*randn(1,Ntap_1)).*sqrt(Power_1/2);
        
        % getting channel impulse response
        h_1(DelayCluster1+1)         =  channel_1;     % cir: channel impulse response
        H_cir_1(:,Sym,ant)           =  h_1;
    end
 end

 for ant = 1:VR2_Antennas
    % rng(0,'twister')
    for Sym = 1:N_realization
        % creating channel in time domain complex entries using number of taps
        channel_2    =  (randn(1,Ntap_2)+ 1i*randn(1,Ntap_2)).*sqrt(Power_2/2);
        
        % getting channel impulse response
        h_2(Time_adjusted_d2+1)      =  channel_2;     % cir: or make this 6
        H_cir_2(:,Sym,ant)           =  h_2;
    end
end

 % figure;waterfall(abs(H_cir_Combined(:,1:100)))
 %    zlabel('Power','FontSize',12,'Fontweight','bold');
 %    ylabel('Time','FontSize',12,'Fontweight','bold');
 %    xlabel('Time','FontSize',12,'Fontweight','bold');
 %    title('Delay clusters combined','FontSize',16,'Fontweight','bold');
 % figure;waterfall(abs(H_cir_1(:,1:100)))
 %    zlabel('Power','FontSize',12,'Fontweight','bold');
 %    ylabel('Time','FontSize',12,'Fontweight','bold');
 %    xlabel('Time','FontSize',12,'Fontweight','bold');
 %    title('Delay cluster 1','FontSize',16,'Fontweight','bold');
 % figure;waterfall(abs(H_cir_2(:,1:100)))
 %    zlabel('Power','FontSize',12,'Fontweight','bold');
 %    ylabel('time','FontSize',12,'Fontweight','bold');
 %    xlabel('Time','FontSize',12,'Fontweight','bold');
 %    title('Delay cluster 2 Time Adjusted','FontSize',16,'Fontweight','bold');
 % figure;waterfall(abs(H_cir_2_D(:,1:100)))
 %    zlabel('Power','FontSize',12,'Fontweight','bold');
 %    ylabel('time','FontSize',12,'Fontweight','bold');
 %    xlabel('Time','FontSize',12,'Fontweight','bold');
 %    title('Delay cluster 2 Normal','FontSize',16,'Fontweight','bold');