%% Demodulation
% BPSK (1), QPSK(2), 8-PSK(3), 16-QAM(4), 64-QAM(6)

function Output_Bits = DeModulator(Modulated_Symbols, Mod_Order)

Output_Bits=zeros(1, length(Modulated_Symbols)*Mod_Order); %Demodulation output

switch Mod_Order
    case 1 % BPSK
        for n= 1 : length(Modulated_Symbols)
            Bit1 = real(Modulated_Symbols(n));
            Output_Bits(n) = Bit1<0;
        end
    
    case 2 % QPSK
        for n= 1 : length(Modulated_Symbols)
            Bit1 = real(Modulated_Symbols(n));
            Bit2 = imag(Modulated_Symbols(n));
            Output_Bits([(Mod_Order*n)-1 (Mod_Order*n)]) = [Bit1<0 Bit2<0];
        end

    case 3 % 8-PSK
        for n= 1 : length(Modulated_Symbols)
            Bit1 = imag(Modulated_Symbols(n));
            Bit2 = real(Modulated_Symbols(n));
            Bit3 = abs(real(Modulated_Symbols(n))) < abs(imag(Modulated_Symbols(n)));
            Output_Bits([(Mod_Order*n)-2 (Mod_Order*n)-1 (Mod_Order*n)]) = [Bit1<0 Bit2<0 Bit3];
        end
        
    case 4 %16-QAM
        for n= 1 : length(Modulated_Symbols)
            Bit1 = real(Modulated_Symbols(n));
            Bit2 = imag(Modulated_Symbols(n));
            Bit3 = (2/sqrt(10)) - abs(real(Modulated_Symbols(n)));
            Bit4 = (2/sqrt(10)) - abs(imag(Modulated_Symbols(n)));
            Output_Bits([(Mod_Order*n)-3 (Mod_Order*n)-2 (Mod_Order*n)-1 (Mod_Order*n)]) = [Bit1<0 Bit2<0 Bit3<0 Bit4<0];
        end
        
    case 6 %64-QAM
        for n= 1 : length(Modulated_Symbols)
            Bit1 = real(Modulated_Symbols(n));
            Bit2 = imag(Modulated_Symbols(n));
            Bit3 = (4/sqrt(42)) - abs(real(Modulated_Symbols(n)));
            Bit4 = (4/sqrt(42)) - abs(imag(Modulated_Symbols(n)));
            if abs(real(Modulated_Symbols(n))) <= (4/sqrt(42))
                Bit5 = abs(real(Modulated_Symbols(n))) - (2/sqrt(42));
            else
                Bit5 = (6/sqrt(42)) - abs(real(Modulated_Symbols(n)));
            end
            if abs(imag(Modulated_Symbols(n))) <= (4/sqrt(42))
                Bit6 = abs(imag(Modulated_Symbols(n))) - (2/sqrt(42));
            else
                Bit6 = (6/sqrt(42)) - abs(imag(Modulated_Symbols(n)));
            end
            Output_Bits([(Mod_Order*n)-5 (Mod_Order*n)-4 (Mod_Order*n)-3 (Mod_Order*n)-2 (Mod_Order*n)-1 (Mod_Order*n)]) = [Bit1<0 Bit2<0 Bit3<0 Bit4<0 Bit5<0 Bit6<0];
        end
end

end