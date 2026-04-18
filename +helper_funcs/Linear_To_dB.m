function value_in_dB = Linear_To_dB(linear_value, flag)
%DB_TO_LINEAR: Converts dB or dBm to linear
%   Detailed explanation goes here
arguments (Input)
    linear_value
    flag
end

arguments (Output)
    value_in_dB
end

if strcmp(flag, 'dB')
    % value_in_dB = 10 * log(linear_value)
    value_in_dB = 10 * log10(linear_value);
elseif strcmp(flag, 'dBm')
    % value_in_dB = 10 * log(linear_value / 1e-3)
    value_in_dB = 10 * log10(linear_value / 1e-3);
else
    error('Invalid flag. Use "dB" or "dBm".');


end