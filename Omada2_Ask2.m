clear; clc; close all;

%% 1. Βασικές παράμετροι συστήματος
M = 4;                  % QPSK
k = log2(M);            % bits ανά σύμβολο
Rs = 5e9;               % symbol rate = 5 Gsym/s
Rb = Rs * k;            % bit rate = 10 Gb/s
alpha = 0.2;            % roll-off factor

Nsym = 5000;            % πλήθος συμβόλων για προσομοίωση

%% 2. Δημιουργία τυχαίων bits
bits = randi([0 1], Nsym * k, 1);

%% 3. Ομαδοποίηση bits ανά 2 για QPSK
bit_pairs = reshape(bits, 2, []).';

%% 4. QPSK Gray mapping
%bit_pairs = reshape(bits, 2, []).';
idx = bi2de(bit_pairs, 'left-msb');

map = [1+1j; -1+1j; 1-1j; -1-1j];
symbols = map(idx + 1) / sqrt(2);
P_avg = mean(abs(symbols).^2);
fprintf('Μέση ισχύς συμβόλων = %.4f\n', P_avg);

%% 5. Κανονικοποίηση ισχύος
%symbols = symbols / sqrt(2);

%% 6. Έλεγχος μέσης ισχύος
P_avg = mean(abs(symbols).^2);

fprintf('Rs = %.2f Gsym/s\n', Rs/1e9);
fprintf('Rb = %.2f Gb/s\n', Rb/1e9);
fprintf('Μέση ισχύς συμβόλων = %.4f\n', P_avg);

%% 7. Διάγραμμα αστερισμού του ιδανικού QPSK
figure;
scatter(real(symbols), imag(symbols), 20, 'filled');
grid on;
axis equal;
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Ιδανικός αστερισμός QPSK');

%% 8. Oversampling και pulse shaping με root-raised-cosine
sps = 8;                 % samples per symbol
span = 10;               % μήκος φίλτρου σε σύμβολα

rrcFilter = rcosdesign(alpha, span, sps, 'sqrt');

% upsampling
symbols_up = upsample(symbols, sps);

% pulse shaping
tx_bb = conv(symbols_up, rrcFilter, 'same');

%% 9. Νέος ρυθμός δειγματοληψίας
Fs = Rs * sps;

fprintf('Fs = %.2f GHz\n', Fs/1e9);

%% 10. Γραφήματα
% (α) Πραγματικό μέρος του baseband σήματος στον χρόνο
figure;
plot(real(tx_bb(1:300)));
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Πραγματικό μέρος του baseband QPSK μετά από RRC pulse shaping');

% (β) Φανταστικό μέρος του baseband σήματος στον χρόνο
figure;
plot(imag(tx_bb(1:300)));
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Φανταστικό μέρος του baseband QPSK μετά από RRC pulse shaping');

% (γ) Αστερισμός στα sample instants
delay = span * sps / 2;
sym_samples = tx_bb(delay+1 : sps : end-delay);

figure;
scatter(real(sym_samples), imag(sym_samples), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός μετά το pulse shaping (στα symbol instants)');

%% 11. Κλιμάκωση του waveform ώστε να αντιστοιχεί σε πραγματική λαμβανόμενη ισχύ
PRX_dBm = -49.1;                         % λαμβανόμενη ισχύς στην είσοδο του δέκτη
PRX_W = 1e-3 * 10^(PRX_dBm/10);         % μετατροπή dBm -> W

P_txbb_norm = mean(abs(tx_bb).^2);      % μέση ισχύς του κανονικοποιημένου waveform

rx_in = tx_bb * sqrt(PRX_W / P_txbb_norm);

P_rx_check_dBm = 10*log10(mean(abs(rx_in).^2)/1e-3);
fprintf('Ισχύς στην είσοδο του δέκτη = %.2f dBm\n', P_rx_check_dBm);

%% 12. 1η βαθμίδα δέκτη: RF BPF ως ισοδύναμη baseband επίδραση
G_BPF_dB = -2;                          % insertion loss 2 dB
G_BPF_lin = 10^(G_BPF_dB/10);

y_bpf = sqrt(G_BPF_lin) * rx_in;

P_bpf_out_dBm = 10*log10(mean(abs(y_bpf).^2)/1e-3);
fprintf('Ισχύς στην έξοδο του BPF = %.2f dBm\n', P_bpf_out_dBm);

%% 13. Διάγραμμα σήματος στην έξοδο του BPF
figure;
plot(real(y_bpf(1:300)));
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Πραγματικό μέρος σήματος στην έξοδο του RF BPF');

%% 14. Αστερισμός στην έξοδο του BPF
sym_bpf = y_bpf(delay+1 : sps : end-delay);

% Κανονικοποίηση μόνο για οπτική σύγκριση του αστερισμού
sym_bpf_norm = sym_bpf / sqrt(mean(abs(sym_bpf).^2));

figure;
scatter(real(sym_bpf_norm), imag(sym_bpf_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του RF BPF');

%% =========================================================
%% 15. 2η βαθμίδα δέκτη: LNA
%% Σε αυτό το βήμα μοντελοποιούμε τον LNA ΜΟΝΟ ως κέρδος.
%% Δεν βάζουμε ακόμα θόρυβο (Noise Figure). Αυτό θα γίνει μετά.
%% =========================================================

G_LNA_dB  = 25;                 % κέρδος LNA από το datasheet
G_LNA_lin = 10^(G_LNA_dB/10);   % μετατροπή από dB σε γραμμική μορφή

% Έξοδος LNA στο complex baseband:
% το πλάτος του σήματος πολλαπλασιάζεται με sqrt(G)
y_lna = sqrt(G_LNA_lin) * y_bpf;

%% Υπολογισμός ισχύος στην έξοδο του LNA
P_lna_out_dBm = 10*log10(mean(abs(y_lna).^2)/1e-3);

fprintf('Ισχύς στην έξοδο του LNA = %.2f dBm\n', P_lna_out_dBm);

%% =========================================================
%% 16. Διάγραμμα σήματος στην έξοδο του LNA
%% Εδώ βλέπουμε το πραγματικό μέρος του waveform μετά την ενίσχυση
%% =========================================================
figure;
plot(real(y_lna(1:300)));
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Πραγματικό μέρος σήματος στην έξοδο του LNA');

%% =========================================================
%% 17. Αστερισμός στην έξοδο του LNA
%% Για οπτική σύγκριση κανονικοποιούμε μόνο το constellation plot
%% ώστε να δούμε αν άλλαξε το σχήμα και όχι μόνο η στάθμη
%% =========================================================
sym_lna = y_lna(delay+1 : sps : end-delay);

sym_lna_norm = sym_lna / sqrt(mean(abs(sym_lna).^2));

figure;
scatter(real(sym_lna_norm), imag(sym_lna_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του LNA');

%% =========================================================
%% 18. Θεωρητικός έλεγχος ισχύος
%% Περιμένουμε:
%% P_LNA,out = P_BPF,out + 25 dB
%% =========================================================
P_lna_expected_dBm = P_bpf_out_dBm + G_LNA_dB;

fprintf('Αναμενόμενη ισχύς στην έξοδο του LNA = %.2f dBm\n', P_lna_expected_dBm);

%% =========================================================
%% 19. Προσθήκη θορύβου LNA μέσω του Noise Figure
%% Εδώ αρχίζουμε να μοντελοποιούμε την υποβάθμιση του SNR
%% =========================================================

NF_LNA_dB = 5;                         % Noise Figure από το datasheet
F_LNA = 10^(NF_LNA_dB/10);             % Noise Factor σε γραμμική μορφή

% Υπολογισμός ισχύος σήματος στην έξοδο του LNA
P_sig_lna = mean(abs(y_lna).^2);

% Προσεγγιστική μοντελοποίηση θορύβου:
% θεωρούμε ότι ο θόρυβος που προστίθεται υποβαθμίζει το SNR
% σύμφωνα με το noise figure του LNA

% Ορίζουμε ένα αρχικό θεωρητικό SNR εισόδου (ενδεικτικά)
SNR_in_dB = 30;                        % ενδεικτική τιμή για αρχικό τεστ

% Το NF μειώνει το SNR στην έξοδο:
SNR_out_LNA_dB = SNR_in_dB - NF_LNA_dB;

% Υπολογισμός ισχύος θορύβου που αντιστοιχεί στο νέο SNR
P_noise_lna = P_sig_lna / (10^(SNR_out_LNA_dB/10));

% Δημιουργία μιγαδικού AWGN
noise_lna = sqrt(P_noise_lna/2) * ...
    (randn(size(y_lna)) + 1j*randn(size(y_lna)));

% Έξοδος LNA με θόρυβο
y_lna_noisy = y_lna + noise_lna;

%% 20. Υπολογισμός SNR στην έξοδο του LNA
P_sig_out = mean(abs(y_lna).^2);
P_noise_out = mean(abs(noise_lna).^2);

SNR_lna_meas_dB = 10*log10(P_sig_out / P_noise_out);

fprintf('SNR στην έξοδο του LNA = %.2f dB\n', SNR_lna_meas_dB);

%% =========================================================
%% 21. Αστερισμός στην έξοδο του LNA με θόρυβο
%% =========================================================
sym_lna_noisy = y_lna_noisy(delay+1 : sps : end-delay);

% Κανονικοποίηση μόνο για οπτική σύγκριση
sym_lna_noisy_norm = sym_lna_noisy / sqrt(mean(abs(sym_lna_noisy).^2));

figure;
scatter(real(sym_lna_noisy_norm), imag(sym_lna_noisy_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του LNA με θόρυβο');

%% =========================================================
%% 22. Υπολογισμός EVM στην έξοδο του LNA
%% =========================================================

% Παίρνουμε τα ιδανικά σύμβολα αναφοράς με ίδιο μήκος
Nvalid = min(length(sym_samples), length(sym_lna_noisy));
ref_sym = sym_samples(1:Nvalid);
rx_sym  = sym_lna_noisy(1:Nvalid);

% Κανονικοποίηση αναφοράς και ληφθέντων συμβόλων
ref_sym = ref_sym / sqrt(mean(abs(ref_sym).^2));
rx_sym  = rx_sym  / sqrt(mean(abs(rx_sym).^2));

evm_rms = sqrt(mean(abs(rx_sym - ref_sym).^2) / mean(abs(ref_sym).^2));
evm_percent = evm_rms * 100;

fprintf('EVM στην έξοδο του LNA = %.2f %%\n', evm_percent);

%% =========================================================
%% 23. 3η βαθμίδα δέκτη: RF Splitter
%% Ο splitter χωρίζει το σήμα σε δύο ίδιους κλάδους (I και Q)
%% και εισάγει απώλεια περίπου -3.5 dB ανά κλάδο
%% =========================================================

G_SPL_dB  = -3.5;                    % συνολική απώλεια ανά κλάδο
G_SPL_lin = 10^(G_SPL_dB/10);

% Έξοδοι splitter στους δύο κλάδους
y_I_split = sqrt(G_SPL_lin) * y_lna_noisy;
y_Q_split = sqrt(G_SPL_lin) * y_lna_noisy;

%% Υπολογισμός ισχύος ανά κλάδο
P_I_split_dBm = 10*log10(mean(abs(y_I_split).^2)/1e-3);
P_Q_split_dBm = 10*log10(mean(abs(y_Q_split).^2)/1e-3);

fprintf('Ισχύς στην έξοδο του RF Splitter (κλάδος I) = %.2f dBm\n', P_I_split_dBm);
fprintf('Ισχύς στην έξοδο του RF Splitter (κλάδος Q) = %.2f dBm\n', P_Q_split_dBm);

%% Θεωρητικός έλεγχος
P_split_expected_dBm = P_lna_out_dBm + G_SPL_dB;
fprintf('Αναμενόμενη ισχύς ανά κλάδο μετά τον splitter = %.2f dBm\n', P_split_expected_dBm);

%% =========================================================
%% 24. Αστερισμός στην έξοδο του splitter (κλάδος I)
%% Οι δύο κλάδοι είναι ίδιοι σε αυτό το στάδιο
%% =========================================================
sym_I_split = y_I_split(delay+1 : sps : end-delay);

sym_I_split_norm = sym_I_split / sqrt(mean(abs(sym_I_split).^2));

figure;
scatter(real(sym_I_split_norm), imag(sym_I_split_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του RF Splitter (κλάδος I)');

%% =========================================================
%% 25. SNR και EVM στην έξοδο του splitter
%% Εφόσον εδώ εφαρμόζουμε μόνο απώλεια, το SNR ιδανικά δεν αλλάζει
%% =========================================================

P_sig_split = mean(abs(sqrt(G_SPL_lin) * y_lna).^2);
P_noise_split = mean(abs(sqrt(G_SPL_lin) * noise_lna).^2);

SNR_split_dB = 10*log10(P_sig_split / P_noise_split);
fprintf('SNR στην έξοδο του RF Splitter = %.2f dB\n', SNR_split_dB);

Nvalid = min(length(sym_samples), length(sym_I_split));
ref_sym = sym_samples(1:Nvalid);
rx_sym  = sym_I_split(1:Nvalid);

ref_sym = ref_sym / sqrt(mean(abs(ref_sym).^2));
rx_sym  = rx_sym  / sqrt(mean(abs(rx_sym).^2));

evm_rms = sqrt(mean(abs(rx_sym - ref_sym).^2) / mean(abs(ref_sym).^2));
evm_percent = evm_rms * 100;

fprintf('EVM στην έξοδο του RF Splitter = %.2f %%\n', evm_percent);
%% =========================================================
%% 26. 4η βαθμίδα δέκτη: Mixer
%% Στο ισοδύναμο complex baseband μοντέλο η ιδανική υποβίβαση
%% συχνότητας θεωρείται ήδη ενσωματωμένη.
%% Άρα εδώ ο mixer μοντελοποιείται μόνο μέσω της conversion loss.
%% =========================================================

G_MIX_dB  = -9;                    % conversion loss από το datasheet
G_MIX_lin = 10^(G_MIX_dB/10);

% Χρησιμοποιούμε έναν από τους δύο ίδιους κλάδους του splitter
% ως είσοδο της βαθμίδας mixer στο ισοδύναμο μοντέλο
y_mix = sqrt(G_MIX_lin) * y_I_split;

%% Υπολογισμός ισχύος στην έξοδο του mixer
P_mix_out_dBm = 10*log10(mean(abs(y_mix).^2)/1e-3);

fprintf('Ισχύς στην έξοδο του Mixer = %.2f dBm\n', P_mix_out_dBm);

%% Θεωρητικός έλεγχος
P_mix_expected_dBm = P_I_split_dBm + G_MIX_dB;
fprintf('Αναμενόμενη ισχύς στην έξοδο του Mixer = %.2f dBm\n', P_mix_expected_dBm);

%% =========================================================
%% 27. Διάγραμμα σήματος στην έξοδο του mixer
%% =========================================================
figure;
plot(real(y_mix(1:300)));
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Πραγματικό μέρος σήματος στην έξοδο του Mixer');

%% =========================================================
%% 28. Αστερισμός στην έξοδο του mixer
%% Κανονικοποιούμε μόνο για οπτική σύγκριση
%% =========================================================
sym_mix = y_mix(delay+1 : sps : end-delay);

sym_mix_norm = sym_mix / sqrt(mean(abs(sym_mix).^2));

figure;
scatter(real(sym_mix_norm), imag(sym_mix_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του Mixer');

%% =========================================================
%% 29. SNR και EVM στην έξοδο του mixer
%% Εφόσον εδώ εφαρμόζουμε μόνο conversion loss,
%% το SNR ιδανικά παραμένει αμετάβλητο
%% =========================================================

% Χωρίζουμε το "καθαρό" σήμα και τον θόρυβο μετά τον splitter
y_I_split_clean = sqrt(G_SPL_lin) * y_lna;
y_I_split_noise = sqrt(G_SPL_lin) * noise_lna;

% Περνάμε και τα δύο από τη conversion loss του mixer
y_mix_clean = sqrt(G_MIX_lin) * y_I_split_clean;
y_mix_noise = sqrt(G_MIX_lin) * y_I_split_noise;

P_sig_mix = mean(abs(y_mix_clean).^2);
P_noise_mix = mean(abs(y_mix_noise).^2);

SNR_mix_dB = 10*log10(P_sig_mix / P_noise_mix);
fprintf('SNR στην έξοδο του Mixer = %.2f dB\n', SNR_mix_dB);

% EVM
Nvalid = min(length(sym_samples), length(sym_mix));
ref_sym = sym_samples(1:Nvalid);
rx_sym  = sym_mix(1:Nvalid);

ref_sym = ref_sym / sqrt(mean(abs(ref_sym).^2));
rx_sym  = rx_sym  / sqrt(mean(abs(rx_sym).^2));

evm_rms = sqrt(mean(abs(rx_sym - ref_sym).^2) / mean(abs(ref_sym).^2));
evm_percent = evm_rms * 100;

fprintf('EVM στην έξοδο του Mixer = %.2f %%\n', evm_percent);

%% =========================================================
%% 30. 5η βαθμίδα δέκτη: LPF
%% Το LPF στο baseband κρατά το χρήσιμο φάσμα μέχρι ~3 GHz
%% και απορρίπτει θόρυβο εκτός ζώνης.
%% Επιπλέον εισάγει απώλεια περίπου -3.1 dB.
%% =========================================================

G_LPF_dB  = -3.1;                   % insertion loss του LPF
G_LPF_lin = 10^(G_LPF_dB/10);

fc_lpf = 3e9;                       % cutoff συχνότητα ~ 3 GHz
Wn = fc_lpf / (Fs/2);               % κανονικοποιημένη συχνότητα για fir1

lpfOrder = 80;                      % τάξη FIR φίλτρου
b_lpf = fir1(lpfOrder, Wn, 'low');  % συντελεστές LPF

% Εφαρμογή LPF στο συνολικό σήμα μετά τον mixer
y_lpf = sqrt(G_LPF_lin) * conv(y_mix, b_lpf, 'same');

%% =========================================================
%% 31. Διαχωρισμός "καθαρού" σήματος και θορύβου μετά το LPF
%% Αυτό χρειάζεται για να υπολογίσουμε SNR στην έξοδο του LPF
%% =========================================================

y_lpf_clean = sqrt(G_LPF_lin) * conv(y_mix_clean, b_lpf, 'same');
y_lpf_noise = sqrt(G_LPF_lin) * conv(y_mix_noise, b_lpf, 'same');

%% =========================================================
%% 32. Ισχύς στην έξοδο του LPF
%% =========================================================
P_lpf_out_dBm = 10*log10(mean(abs(y_lpf).^2)/1e-3);
fprintf('Ισχύς στην έξοδο του LPF = %.2f dBm\n', P_lpf_out_dBm);

P_lpf_expected_dBm = P_mix_out_dBm + G_LPF_dB;
fprintf('Αναμενόμενη ισχύς στην έξοδο του LPF = %.2f dBm\n', P_lpf_expected_dBm);

%% =========================================================
%% 33. Διάγραμμα σήματος στην έξοδο του LPF
%% =========================================================
figure;
plot(real(y_lpf(1:300)));
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Πραγματικό μέρος σήματος στην έξοδο του LPF');

%% =========================================================
%% 34. Αστερισμός στην έξοδο του LPF
%% Κανονικοποίηση μόνο για οπτική σύγκριση
%% =========================================================
sym_lpf = y_lpf(delay+1 : sps : end-delay);

sym_lpf_norm = sym_lpf / sqrt(mean(abs(sym_lpf).^2));

figure;
scatter(real(sym_lpf_norm), imag(sym_lpf_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του LPF');

%% =========================================================
%% 35. Υπολογισμός SNR στην έξοδο του LPF
%% =========================================================
P_sig_lpf = mean(abs(y_lpf_clean).^2);
P_noise_lpf = mean(abs(y_lpf_noise).^2);

SNR_lpf_dB = 10*log10(P_sig_lpf / P_noise_lpf);
fprintf('SNR στην έξοδο του LPF = %.2f dB\n', SNR_lpf_dB);

%% =========================================================
%% 36. Υπολογισμός EVM στην έξοδο του LPF
%% =========================================================
Nvalid = min(length(sym_samples), length(sym_lpf));
ref_sym = sym_samples(1:Nvalid);
rx_sym  = sym_lpf(1:Nvalid);

ref_sym = ref_sym / sqrt(mean(abs(ref_sym).^2));
rx_sym  = rx_sym  / sqrt(mean(abs(rx_sym).^2));

evm_rms = sqrt(mean(abs(rx_sym - ref_sym).^2) / mean(abs(ref_sym).^2));
evm_percent = evm_rms * 100;

fprintf('EVM στην έξοδο του LPF = %.2f %%\n', evm_percent);

%% =========================================================
%% 37. 6η βαθμίδα δέκτη: Baseband VGA
%% Το VGA αυξάνει τη στάθμη του baseband σήματος πριν από τον ADC.
%% Σε αυτό το βήμα βάζουμε μόνο το gain του VGA.
%% =========================================================

G_VGA_dB  = 16;                    % gain του VGA από το operating point
G_VGA_lin = 10^(G_VGA_dB/10);

% Έξοδος VGA
y_vga = sqrt(G_VGA_lin) * y_lpf;

%% =========================================================
%% 38. Ισχύς στην έξοδο του VGA
%% =========================================================
P_vga_out_dBm = 10*log10(mean(abs(y_vga).^2)/1e-3);
fprintf('Ισχύς στην έξοδο του VGA = %.2f dBm\n', P_vga_out_dBm);

P_vga_expected_dBm = P_lpf_out_dBm + G_VGA_dB;
fprintf('Αναμενόμενη ισχύς στην έξοδο του VGA = %.2f dBm\n', P_vga_expected_dBm);

%% =========================================================
%% 39. Διάγραμμα σήματος στην έξοδο του VGA
%% =========================================================
figure;
plot(real(y_vga(1:300)));
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Πραγματικό μέρος σήματος στην έξοδο του VGA');

%% =========================================================
%% 40. Αστερισμός στην έξοδο του VGA
%% Κανονικοποίηση μόνο για οπτική σύγκριση
%% =========================================================
sym_vga = y_vga(delay+1 : sps : end-delay);

sym_vga_norm = sym_vga / sqrt(mean(abs(sym_vga).^2));

figure;
scatter(real(sym_vga_norm), imag(sym_vga_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του VGA');

%% =========================================================
%% 41. SNR και EVM στην έξοδο του VGA
%% Εφόσον εδώ εφαρμόζουμε μόνο gain, το SNR ιδανικά δεν αλλάζει
%% =========================================================

y_vga_clean = sqrt(G_VGA_lin) * y_lpf_clean;
y_vga_noise = sqrt(G_VGA_lin) * y_lpf_noise;

P_sig_vga = mean(abs(y_vga_clean).^2);
P_noise_vga = mean(abs(y_vga_noise).^2);

SNR_vga_dB = 10*log10(P_sig_vga / P_noise_vga);
fprintf('SNR στην έξοδο του VGA = %.2f dB\n', SNR_vga_dB);

Nvalid = min(length(sym_samples), length(sym_vga));
ref_sym = sym_samples(1:Nvalid);
rx_sym  = sym_vga(1:Nvalid);

ref_sym = ref_sym / sqrt(mean(abs(ref_sym).^2));
rx_sym  = rx_sym  / sqrt(mean(abs(rx_sym).^2));

evm_rms = sqrt(mean(abs(rx_sym - ref_sym).^2) / mean(abs(ref_sym).^2));
evm_percent = evm_rms * 100;

fprintf('EVM στην έξοδο του VGA = %.2f %%\n', evm_percent);

%% =========================================================
%% 42. 7η βαθμίδα δέκτη: ADC
%% Ο ADC είναι dual, άρα στο complex baseband μοντέλο
%% κβαντίζουμε ξεχωριστά το I και το Q κανάλι.
%% =========================================================

ADC_bits = 12;                         % ανάλυση ADC
P_FS_ADC_dBm = 3.9;                   % full-scale input power από το datasheet
P_FS_ADC_W = 1e-3 * 10^(P_FS_ADC_dBm/10);

% Για κάθε πραγματικό κανάλι (I ή Q) θεωρούμε full-scale ημιτονοειδές
A_FS_rms  = sqrt(P_FS_ADC_W);         % RMS πλάτος που αντιστοιχεί στο full-scale power
A_FS_peak = sqrt(2) * A_FS_rms;       % peak amplitude για ημίτονο

% Βήμα κβαντοποίησης
Delta = 2*A_FS_peak / (2^ADC_bits - 1);

%% =========================================================
%% 43. Διαχωρισμός σε I και Q πριν τον ADC
%% =========================================================
I_adc_in = real(y_vga);
Q_adc_in = imag(y_vga);

%% =========================================================
%% 44. Clipping στα όρια full-scale
%% Αν το σήμα ξεπερνά τα όρια του ADC, κόβεται
%% =========================================================
I_clip = min(max(I_adc_in, -A_FS_peak), A_FS_peak);
Q_clip = min(max(Q_adc_in, -A_FS_peak), A_FS_peak);

%% =========================================================
%% 45. Κβαντοποίηση 12-bit
%% =========================================================
I_adc = Delta * round(I_clip / Delta);
Q_adc = Delta * round(Q_clip / Delta);

% Επανασύνθεση μιγαδικού baseband σήματος
y_adc = I_adc + 1j*Q_adc;

%% =========================================================
%% 46. Έλεγχος ισχύος και στάθμης στον ADC
%% =========================================================
P_adc_out_dBm = 10*log10(mean(abs(y_adc).^2)/1e-3);
fprintf('Ισχύς στην έξοδο του ADC (μετά την κβαντοποίηση) = %.2f dBm\n', P_adc_out_dBm);

P_I_adc_in_dBm = 10*log10(mean(I_adc_in.^2)/1e-3);
P_Q_adc_in_dBm = 10*log10(mean(Q_adc_in.^2)/1e-3);

fprintf('Ισχύς στο κανάλι I πριν τον ADC = %.2f dBm\n', P_I_adc_in_dBm);
fprintf('Ισχύς στο κανάλι Q πριν τον ADC = %.2f dBm\n', P_Q_adc_in_dBm);

fprintf('Full-scale input power ADC = %.2f dBm\n', P_FS_ADC_dBm);

Margin_I_dB = P_FS_ADC_dBm - P_I_adc_in_dBm;
Margin_Q_dB = P_FS_ADC_dBm - P_Q_adc_in_dBm;

fprintf('Περιθώριο από full-scale στο κανάλι I = %.2f dB\n', Margin_I_dB);
fprintf('Περιθώριο από full-scale στο κανάλι Q = %.2f dB\n', Margin_Q_dB);

%% =========================================================
%% 47. Διάγραμμα χρόνου: είσοδος και έξοδος ADC (κανάλι I)
%% =========================================================
Nplot = 150;

figure;
plot(I_adc_in(1:Nplot), 'LineWidth', 1.1);
hold on;
stairs(I_adc(1:Nplot), 'LineWidth', 1.1);
grid on;
xlabel('Δείγματα');
ylabel('Πλάτος');
title('Κανάλι I πριν και μετά τον ADC');
legend('Είσοδος ADC', 'Έξοδος ADC (κβαντισμένη)');

%% =========================================================
%% 48. Αστερισμός στην έξοδο του ADC
%% =========================================================
sym_adc = y_adc(delay+1 : sps : end-delay);

sym_adc_norm = sym_adc / sqrt(mean(abs(sym_adc).^2));

figure;
scatter(real(sym_adc_norm), imag(sym_adc_norm), 20, 'filled');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Αστερισμός στην έξοδο του ADC');

%% =========================================================
%% 49. Effective SNR στην έξοδο του ADC
%% Το σφάλμα τώρα περιλαμβάνει:
%% - αναλογικό θόρυβο προηγούμενων βαθμίδων
%% - κβαντοποίηση
%% - πιθανό clipping
%% =========================================================

P_sig_adc = mean(abs(y_vga_clean).^2);
P_err_adc = mean(abs(y_adc - y_vga_clean).^2);

SNR_adc_eff_dB = 10*log10(P_sig_adc / P_err_adc);
fprintf('Effective SNR στην έξοδο του ADC = %.2f dB\n', SNR_adc_eff_dB);

%% =========================================================
%% 50. EVM στην έξοδο του ADC
%% =========================================================
Nvalid = min(length(sym_samples), length(sym_adc));
ref_sym = sym_samples(1:Nvalid);
rx_sym  = sym_adc(1:Nvalid);

ref_sym = ref_sym / sqrt(mean(abs(ref_sym).^2));
rx_sym  = rx_sym  / sqrt(mean(abs(rx_sym).^2));

evm_rms = sqrt(mean(abs(rx_sym - ref_sym).^2) / mean(abs(ref_sym).^2));
evm_percent = evm_rms * 100;

fprintf('EVM στην έξοδο του ADC = %.2f %%\n', evm_percent);