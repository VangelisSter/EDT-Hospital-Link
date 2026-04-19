function evm_percent = Calculate_EVM(y_clean,noise)
%CALCULATE_EVM Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    y_clean
    noise
end

arguments (Output)
    evm_percent
end

rms_noise = sqrt(mean(abs(noise).^2));
rms_signal = sqrt(mean(abs(y_clean).^2));

evm_percent = (rms_noise / rms_signal) * 100;
end