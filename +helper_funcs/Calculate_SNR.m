function SNR = Calculate_SNR(y,y_noise)
%CALCULATE_SNR Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    y
    y_noise
end

arguments (Output)
    SNR
end

P_sig_out = mean(abs(y).^2);
P_noise_out = mean(abs(y_noise).^2);

SNR = 10*log10(P_sig_out / P_noise_out);
end