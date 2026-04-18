function linear_value = dB_To_Linear(value_in_dB, flag)
%DB_TO_LINEAR: Converts dB or dBm to linear
%   Detailed explanation goes here
arguments (Input)
    value_in_dB
    flag
end

arguments (Output)
    linear_value
end

if strcmp(flag, 'dB')
    % value_in_dB = 10 * log(linear_value)
    linear_value = 10^(value_in_dB / 10);
elseif strcmp(flag, 'dBm')
    % value_in_dB = 10 * log(linear_value / 1e-3)
    linear_value = 1e-3 * 10^(value_in_dB  / 10);
else
    error('Invalid flag. Use "dB" or "dBm".');


end