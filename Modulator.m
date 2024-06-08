%% Modulation
% BPSK (1), QPSK(2), 8-PSK(3), 16-QAM(4), 64-QAM(6)

function Modulated_Symbols = Modulator(Input_Bits, Mod_Order)

Modulated_Symbols = zeros(1, length(Input_Bits)/Mod_Order); %Modulation output

switch Mod_Order
    
    case 1  % BPSK
        Modulated_Symbols = -1*(2*(Input_Bits)-1);
        
    case 2 % QPSK
        for n= 1 : Mod_Order : length(Input_Bits)
            I= 1; Q= 1;
            if Input_Bits(n) == 1 %Check for real part sign
                I= -I;
            end
            if Input_Bits(n+1) == 1 %Check for imaginary part sign
                Q= -Q;
            end
            Modulated_Symbols(floor(n/Mod_Order)+1) = (I+1i*Q)/sqrt(2);
        end
        
    case 3 % 8-PSK
        for n= 1 : Mod_Order : length(Input_Bits)
            I= 1; Q= 1;
            if Input_Bits(n) == 1    %Check for imaginary part sign
                Q= -Q;
            end
            if Input_Bits(n+1) == 1  %Check for real part sign
                I= -I;
            end
            if Input_Bits(n+2) == 1  %Check for phase
                Q= 3*Q;
            else
                I= 3*I;
            end
            Modulated_Symbols(floor(n/Mod_Order)+1) = (I+1i*Q)/sqrt(10);
        end
        
    case 4 % 16-QAM
        for n= 1 : Mod_Order : length(Input_Bits)
            I= 1;Q= 1;
            if Input_Bits(n) == 1    %Check for real part sign
                I= -I;
            end
            if Input_Bits(n+1) == 1  %Check for imaginary part sign
                Q= -Q;
            end
            if Input_Bits(n+2) == 1  %Check for real part value
                I= 3*I;
            end
            if Input_Bits(n+3) == 1  %Check for imaginary part value
                Q= 3*Q;
            end
            Modulated_Symbols(floor(n/Mod_Order)+1) = (I+1i*Q)/sqrt(10); 
        end
        
    case 6 %64-QAM
        for n = 1 : Mod_Order : length(Input_Bits)
            I= 1; Q= 1;
            if Input_Bits(n) == 1    %Check for real part sign
                I= -I;
            end
            if Input_Bits(n+1) == 1  %Check for imaginary part sign
                Q= -Q;
            end
            if Input_Bits(n+2) == 0  %Check for real part value
                if Input_Bits(n+4) == 0
                    I= 3*I;
                end
            else
                if Input_Bits(n+4) == 0
                    I= 5*I;
                else
                    I= 7*I;
                end
            end
            if Input_Bits(n+3) == 0  %Check for imaginary part value
                if Input_Bits(n+5) == 0
                    Q= 3*Q;
                end
            else
                if Input_Bits(n+5) == 0
                    Q= 5*Q;
                else
                    Q= 7*Q;
                end
            end
            Modulated_Symbols(floor(n/Mod_Order)+1) = (I+1i*Q)/sqrt(42);
        end
end

end