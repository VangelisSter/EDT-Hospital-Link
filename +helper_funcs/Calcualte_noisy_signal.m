function [noisy_signal,SNR_out_dB] = Calcualte_noisy_signal(NF_dB,signal_in, SNR_in_dB)
%CALCUALTE_NOISY_SIGNAL Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    NF_dB
    signal_in
    SNR_in_dB
end

arguments (Output)
    noisy_signal
    SNR_out_dB
end

% Υπολογισμός ισχύος σήματος στην έξοδο του stage
P_sig_in = mean(abs(signal_in).^2);

% Προσεγγιστική μοντελοποίηση θορύβου:
% θεωρούμε ότι ο θόρυβος που προστίθεται υποβαθμίζει το SNR

% Το NF μειώνει το SNR στην έξοδο:
SNR_out_dB = SNR_in_dB - NF_dB;

% Υπολογισμός ισχύος θορύβου που αντιστοιχεί στο νέο SNR
P_noise = P_sig_in / (10^(SNR_out_dB/10));

% Δημιουργία μιγαδικού AWGN
noise = sqrt(P_noise/2) * ...
    (randn(size(signal_in)) + 1j*randn(size(signal_in)));

% Έξοδος LNA με θόρυβο
noisy_signal = signal_in + noise;



end