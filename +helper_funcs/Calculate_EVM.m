function evm_percent = Calculate_EVM(sym_samples,sym_noisy)
%CALCULATE_EVM Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    sym_samples
    sym_noisy
end

arguments (Output)
    evm_percent
end

% Παίρνουμε τα ιδανικά σύμβολα αναφοράς με ίδιο μήκος
Nvalid = min(length(sym_samples), length(sym_noisy));
ref_sym = sym_samples(1:Nvalid);
rx_sym  = sym_noisy(1:Nvalid);

% Κανονικοποίηση αναφοράς και ληφθέντων συμβόλων
ref_sym = ref_sym / sqrt(mean(abs(ref_sym).^2));
rx_sym  = rx_sym  / sqrt(mean(abs(rx_sym).^2));

evm_rms = sqrt(mean(abs(rx_sym - ref_sym).^2) / mean(abs(ref_sym).^2));
evm_percent = evm_rms * 100;
end