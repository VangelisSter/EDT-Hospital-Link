function [noisy_signal, noise] = Calcualte_noisy_signal(signal_in, SNR_in_dB, SNR_out_dB, flag)
%CALCUALTE_NOISY_SIGNAL Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    signal_in
    SNR_in_dB
    SNR_out_dB
    flag
end

arguments (Output)
    noisy_signal 
    noise
end

% Υπολογισμός ισχύος σήματος στην έξοδο του stage
P_sig_out = mean(abs(signal_in).^2);

% Προσεγγιστική μοντελοποίηση θορύβου:
% θεωρούμε ότι ο θόρυβος που προστίθεται υποβαθμίζει το SNR

P_noise_existing = P_sig_out / (10^(SNR_in_dB/10));

% Υπολογισμός ισχύος θορύβου που αντιστοιχεί στο νέο SNR
P_noise_total = P_sig_out / (10^(SNR_out_dB/10));

P_noise_to_add = P_noise_total - P_noise_existing;

P_noise_to_add = max(P_noise_to_add, 0);


% Δημιουργία μιγαδικού AWGN
if strcmp(flag,'complex')
    noise = sqrt(P_noise_to_add/2) * ...
        (randn(size(signal_in)) + 1j*randn(size(signal_in)));
else
    noise = sqrt(P_noise_to_add/2) * ...
        (randn(size(signal_in)));
end

% Έξοδος LNA με θόρυβο
noisy_signal = signal_in + noise;



end