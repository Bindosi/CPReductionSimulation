function [H_cir_1, Lch_1, H_cir_2, Lch_2, Lch_Combined] = Maximal_Ratio_Combining(PowerdBCluser1,PowerdBCluser2)
Power1  =  [-3 -9 -19 -25 -31];     % Cluster 1 tap power profile 'dB'
Power2  =  [-1 -5 -11 -19 -26 33];  % Cluster 1 tap power profile 'dB'

VR1_antennas      =     64;
VR2_antennas      =     40;

Ntap1   =  length(Power1);
Ntap2   =  length(Power2);

channel1            =  (randn(1,Ntap1)  + 1i*randn(1,Ntap1))  .*sqrt(Power1/2);
channel2            =  (randn(1,Ntap2)  + 1i*randn(1,Ntap2))  .*sqrt(Power2/2);

Channel_Vector      = [channel1,channel2];
Norm_Channel_Vector = norm(Channel_Vector);
Combining_Vector    = conj(Channel_Vector'/Norm_Channel_Vector);

