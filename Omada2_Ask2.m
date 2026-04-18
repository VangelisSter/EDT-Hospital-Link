clear;
clc;
close all
%Parameters
M = 4;
k = log2(M);
sps = 32;
fc = 70e9;
BW = 6e9;
rolloff = 0.2;
filtlen = 20; % Filter cut-off at 20 syms
Rs = BW/(1 + rolloff);
Rb = Rs * k;            % bit rate = 10 Gb/s
Fs = sps * Rs; %συχνοτητα δειγματοληψιας, πρεπει τουλαχιστον 2 επι το symbol rate
N_syms = 10000;
bits = randi([0,1], 1 , N_syms * k);
idx = bi2de(reshape(bits, k , []).' , 'left-msb');
% Rotate PSK by pi/4 so the values are not on the axis
syms = pskmod(idx, M, pi/4);
% filtlen*sps + 1 = 81 so 81 coefficients, odd because rc is
% symmetrical around 0
rrc = rcosdesign( rolloff, filtlen, sps, 'sqrt');
rrc = rrc / sqrt(sum(rrc.^2));
% umpsample so each symbol becomes sps = 8 samples
% rrc in each symbol
tx_shaped = upfirdn(syms, rrc, sps);
delay = (filtlen * sps) / 2;
% Keep only after the transient delay introduced by the filter
tx_shaped_valid = tx_shaped(delay + 1: end - delay);
tx_shaped_valid = tx_shaped_valid / sqrt(mean(abs(tx_shaped_valid).^2));
% real_tx = real(tx_shaped_valid) / max(real(tx_shaped_valid)) * max(real(syms));
% imag_tx = imag(tx_shaped_valid) / max(imag(tx_shaped_valid)) * max(imag(syms));
% tx_shaped_valid = real_tx + 1i .* imag_tx;
% Take the I and Q of the complex values
tx_I = real(tx_shaped_valid');
tx_Q = imag(tx_shaped_valid');
N = length(tx_shaped_valid); % total number of samples
t = (0:N-1)/Fs; % Time between consecutive samples

%% 6. Έλεγχος μέσης ισχύος
P_avg = mean(abs(syms).^2);

fprintf('Rs = %.2f Gsym/s\n', Rs/1e9);
fprintf('Rb = %.2f Gb/s\n', Rb/1e9);
fprintf('Μέση ισχύς συμβόλων = %.4f\n', P_avg);

%% 7. Διάγραμμα αστερισμού του ιδανικού QPSK
helper_funcs.Plot_Constellation(syms, "Ιδανικός αστερισμός QPSK");

%% 9. Νέος ρυθμός δειγματοληψίας
Fs = Rs * sps;

fprintf('Fs = %.2f GHz\n', Fs/1e9);

%% 10. Γραφήματα
% (α) Πραγματικό μέρος του baseband σήματος στον χρόνο
helper_funcs.Plot_RF_Signal(t, tx_I, "Πραγματικό μέρος του baseband QPSK μετά από RRC pulse shaping");

% (β) Φανταστικό μέρος του baseband σήματος στον χρόνο
helper_funcs.Plot_RF_Signal(t, tx_Q, "Φανταστικό μέρος του baseband QPSK μετά από RRC pulse shaping");

% (γ) Αστερισμός στα sample instants
%delay = span * sps / 2;
sym_samples = tx_shaped_valid(delay+1 : sps : end-delay);
helper_funcs.Plot_Constellation(sym_samples,"Αστερισμός μετά το pulse shaping (στα symbol instants)");

%% 11. Κλιμάκωση του waveform ώστε να αντιστοιχεί σε πραγματική λαμβανόμενη ισχύ
SNRin_dB = 20;% SNR at the receiver
k = 1.38e-23;  % Boltzman constant
T_A = 150;     % In Kelvin
Nin = k * T_A * BW; % The formula for the noise
Nin_dBm = 10 * log10(Nin / 1e-3);

% We know that SNR = 10log(S/N) so S = (10 ^ (SNR / 10)) * N
Sin = (10 ^ (SNRin_dB / 10)) * Nin;
SNR_in_test = 10 * log10(Sin / Nin);
PRX_dBm = 10 * log10(Sin / 1e-3);
PRX_W = 1e-3 * 10^(PRX_dBm/10);         % Convert dBm -> W

P_txbb_norm = mean(abs(tx_shaped_valid).^2);      % Mean Power of normalized waveform

rx_in = tx_shaped_valid * sqrt(PRX_W * 50 / P_txbb_norm);

P_rx_check_dBm = 10*log10(mean(abs(rx_in).^2)/1e-3/50);
fprintf('Ισχύς στην είσοδο του δέκτη = %.2f dBm\n', P_rx_check_dBm);

%% 12. 1η βαθμίδα δέκτη: RF BPF ως ισοδύναμη baseband επίδραση
G_BPF_dB = -1.5;                          % insertion loss 2 dB
G_BPF_lin = 10^(G_BPF_dB/10);
NF_BPF_dB = - G_BPF_dB; % Not mentioned in the datasheet
NF_BPF = helper_funcs.dB_To_Linear(NF_BPF_dB, 'dB');

y_bpf = sqrt(G_BPF_lin) * rx_in;
[y_bpf_noisy, SNR_out_BPF_dB, bpf_noise] = helper_funcs.Calcualte_noisy_signal(NF_BPF_dB, y_bpf, SNRin_dB);

P_bpf_out_dBm = 10*log10(mean(abs(y_bpf).^2)/1e-3/50);
P_bpf_expected_dBm = PRX_dBm + G_BPF_dB;
fprintf('Expected Power after BPF: %.2f dBm\n', P_bpf_expected_dBm);
fprintf('Power after BPF: %.2f dBm\n', P_bpf_out_dBm);

% Real Part of the signal after the BPF with noise
helper_funcs.Plot_RF_Signal(t, real(y_bpf_noisy), "Πραγματικό μέρος σήματος στην έξοδο του RF BPF");

% Constellation after the BPF with noise
sym_bpf = y_bpf_noisy(delay+1 : sps : end-delay);

% We normalize so that we can properly compare
sym_bpf_norm = sym_bpf / sqrt(mean(abs(sym_bpf).^2));

helper_funcs.Plot_Constellation(sym_bpf_norm, "Αστερισμός στην έξοδο του RF BPF");

SNR_out_BPF_calc_dB = helper_funcs.Calculate_SNR(y_bpf, bpf_noise);
fprintf('SNR στην έξοδο του BPF = %.2f dB\n', SNR_out_BPF_calc_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_lna);
fprintf('EVM στην έξοδο του LNA = %.2f %%\n', evm_percent);

%% =========================================================
% 15. 2η βαθμίδα δέκτη: LNA
% =========================================================

G_LNA_dB  = 25;               
G_LNA_lin = 10^(G_LNA_dB/10);
NF_LNA_dB = 5;
NF_LNA = helper_funcs.dB_To_Linear(NF_LNA_dB, 'dB');

% Έξοδος LNA στο complex baseband:
% το πλάτος του σήματος πολλαπλασιάζεται με sqrt(G)
y_lna = sqrt(G_LNA_lin) * y_bpf;
[y_lna_noisy, SNR_out_LNA] = helper_funcs.Calcualte_noisy_signal(NF_LNA_dB, y_lna, SNR_out_BPF_dB);

% Real Part of the signal after the LNA with noise
helper_funcs.Plot_RF_Signal(t, real(y_lna_noisy), "Real part of signal after LNA");

% Constellation after the LNA with noise
sym_lna = y_lna_noisy(delay+1 : sps : end-delay);

sym_lna_noisy_norm = sym_lna / sqrt(mean(abs(sym_lna).^2));

helper_funcs.Plot_Constellation(sym_lna_noisy_norm, "Αστερισμός στην έξοδο του LNA");

% Power calculations after LNA
% We expect P_LNA,out = P_BPF,out + 25 dB
P_lna_out_dBm = 10*log10(mean(abs(y_lna).^2)/1e-3/50);
P_lna_expected_dBm = P_bpf_out_dBm + G_LNA_dB;

fprintf('Expected Power after LNA: %.2f dBm\n', P_lna_expected_dBm);
fprintf('Power after LNA: %.2f dBm\n', P_lna_out_dBm);

% Now we calculate the SNR and EVM after LNA

SNR_lna_meas_dB = helper_funcs.Calculate_SNR(y_lna, y_lna_noisy);
fprintf('SNR στην έξοδο του LNA = %.2f dB\n', SNR_lna_meas_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_lna);
fprintf('EVM στην έξοδο του LNA = %.2f %%\n', evm_percent);

%% =========================================================
% 23. 3η βαθμίδα δέκτη: RF Splitter
% Ο splitter χωρίζει το σήμα σε δύο ίδιους κλάδους (I και Q)
% και εισάγει απώλεια περίπου -3.5 dB ανά κλάδο
% =========================================================

G_PS_dB  = -3.5;                    % συνολική απώλεια ανά κλάδο
G_PS_lin = 10^(G_PS_dB/10);
NF_PS_dB = - G_PS_dB; % Not mentioned in the datasheet
NF_PS = helper_funcs.dB_To_Linear(NF_PS_dB, 'dB');

% Έξοδοι splitter στους δύο κλάδους
y_I_split = sqrt(G_PS_lin) * y_lna;
y_Q_split = sqrt(G_PS_lin) * y_lna;
[y_I_split_noisy, SNR_out_I_split] = helper_funcs.Calcualte_noisy_signal(NF_PS_dB, y_lna_noisy, SNR_out_LNA);
[y_Q_split_noisy, SNR_out_Q_split] = helper_funcs.Calcualte_noisy_signal(NF_PS_dB, y_lna_noisy, SNR_out_LNA);

helper_funcs.Plot_RF_Signal(t, real(y_I_split_noisy), "Real part of signal after Power Splitter (Branch I)");

sym_I_split = y_I_split_noisy(delay+1 : sps : end-delay);

sym_I_split_norm = sym_I_split / sqrt(mean(abs(sym_I_split).^2));

helper_funcs.Plot_Constellation(sym_I_split_norm, "Constellation after RF Splitter (branch I)");

% Υπολογισμός ισχύος ανά κλάδο
P_I_split_dBm = 10*log10(mean(abs(y_I_split).^2)/1e-3);
P_Q_split_dBm = 10*log10(mean(abs(y_Q_split).^2)/1e-3);
P_split_expected_dBm = P_lna_out_dBm + G_PS_dB;

fprintf('Expected Power per branch after Power Splitter: %.2f dBm\n', P_split_expected_dBm);
fprintf('Power at branch I after Power Splitter: %.2f dBm\n', P_I_split_dBm);
fprintf('Power at branch Q after Power Splitter: %.2f dBm\n', P_Q_split_dBm);

% Now we calculate the SNR and EVM after Power Spliiter

SNR_PS_dB = helper_funcs.Calculate_SNR(y_I_split, y_I_split_noisy);
fprintf('SNR στην έξοδο του RF Splitter = %.2f dB\n', SNR_PS_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_I_split);
fprintf('EVM στην έξοδο του RF Splitter = %.2f %%\n', evm_percent);

%% =========================================================
% 26. 4η βαθμίδα δέκτη: Mixer
% Στο ισοδύναμο complex baseband μοντέλο η ιδανική υποβίβαση
% συχνότητας θεωρείται ήδη ενσωματωμένη.
% Άρα εδώ ο mixer μοντελοποιείται μόνο μέσω της conversion loss.
% =========================================================

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
% 27. Διάγραμμα σήματος στην έξοδο του mixer
% =========================================================

helper_funcs.Plot_RF_Signal(t, real(y_mix), "Πραγματικό μέρος σήματος στην έξοδο του Mixer");

%% =========================================================
% 28. Αστερισμός στην έξοδο του mixer
% Κανονικοποιούμε μόνο για οπτική σύγκριση
%%=========================================================
sym_mix = y_mix(delay+1 : sps : end-delay);

sym_mix_norm = sym_mix / sqrt(mean(abs(sym_mix).^2));

helper_funcs.Plot_Constellation(sym_mix_norm, "Αστερισμός στην έξοδο του Mixer");

%% =========================================================
% 29. SNR και EVM στην έξοδο του mixer
% Εφόσον εδώ εφαρμόζουμε μόνο conversion loss,
% το SNR ιδανικά παραμένει αμετάβλητο
% =========================================================

% Χωρίζουμε το "καθαρό" σήμα και τον θόρυβο μετά τον splitter
y_I_split_clean = sqrt(G_PS_lin) * y_lna;
y_I_split_noise = sqrt(G_PS_lin) * noise_lna;

% Περνάμε και τα δύο από τη conversion loss του mixer
y_mix_clean = sqrt(G_MIX_lin) * y_I_split_clean;
y_mix_noise = sqrt(G_MIX_lin) * y_I_split_noise;

SNR_mix_dB = helper_funcs.Calculate_SNR(y_mix_clean, y_mix_noise);
fprintf('SNR στην έξοδο του Mixer = %.2f dB\n', SNR_mix_dB);

% EVM
evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_mix);
fprintf('EVM στην έξοδο του Mixer = %.2f %%\n', evm_percent);

%% =========================================================
% 30. 5η βαθμίδα δέκτη: LPF
% Το LPF στο baseband κρατά το χρήσιμο φάσμα μέχρι ~3 GHz
% και απορρίπτει θόρυβο εκτός ζώνης.
% Επιπλέον εισάγει απώλεια περίπου -3.1 dB.
% =========================================================

G_LPF_dB  = -3.1;                   % insertion loss του LPF
G_LPF_lin = 10^(G_LPF_dB/10);

fc_lpf = 3e9;                       % cutoff συχνότητα ~ 3 GHz
Wn = fc_lpf / (Fs/2);               % κανονικοποιημένη συχνότητα για fir1

lpfOrder = 80;                      % τάξη FIR φίλτρου
b_lpf = fir1(lpfOrder, Wn, 'low');  % συντελεστές LPF

% Εφαρμογή LPF στο συνολικό σήμα μετά τον mixer
y_lpf = sqrt(G_LPF_lin) * conv(y_mix, b_lpf, 'same');
y_lpf_rrc = conv(y_lpf, rrc, 'same');
power_before = mean(abs(y_lpf).^2);
power_after = mean(abs(y_lpf_rrc).^2);

% Optional: Re-scale if there is a discrepancy
y_lpf_rrc = y_lpf_rrc * sqrt(power_before / power_after);
helper_funcs.Plot_RF_Signal(t, real(y_lpf_rrc), "LPF Signal with RRC before normalization");
%y_lpf_rrc = y_lpf_rrc / sqrt(mean(abs(y_lpf_rrc).^2));
%%
% Effect of RRC
figure;
[pxx, f] = pwelch(y_lpf, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Φάσμα Baseband LPF χωρίς RRC');
grid on;
xlim([-BW/(2*1e9) - 2  BW/(2*1e9) + 2]);

figure;
[pxx, f] = pwelch(rrc, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Φάσμα Baseband RRC');
grid on;
xlim([-BW/(2*1e9) - 2  BW/(2*1e9) + 2]);

figure;
[pxx, f] = pwelch(y_lpf_rrc, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Φάσμα Baseband LPF me RRC');
grid on;
xlim([-BW/(2*1e9) - 2  BW/(2*1e9) + 2]);
%% =========================================================
% 31. Διαχωρισμός "καθαρού" σήματος και θορύβου μετά το LPF
% Αυτό χρειάζεται για να υπολογίσουμε SNR στην έξοδο του LPF
% =========================================================

y_lpf_clean = sqrt(G_LPF_lin) * conv(y_mix_clean, b_lpf, 'same');
y_lpf_noise = sqrt(G_LPF_lin) * conv(y_mix_noise, b_lpf, 'same');

%% =========================================================
% 32. Ισχύς στην έξοδο του LPF
% =========================================================
P_lpf_out_dBm = 10*log10(mean(abs(y_lpf_rrc).^2)/1e-3);
fprintf('Ισχύς στην έξοδο του LPF = %.2f dBm\n', P_lpf_out_dBm);

P_lpf_expected_dBm = P_mix_out_dBm + G_LPF_dB;
fprintf('Αναμενόμενη ισχύς στην έξοδο του LPF = %.2f dBm\n', P_lpf_expected_dBm);

%% =========================================================
% 33. Διάγραμμα σήματος στην έξοδο του LPF
% =========================================================

helper_funcs.Plot_RF_Signal(t, real(y_lpf), "Πραγματικό μέρος σήματος στην έξοδο του LPF χωρίς RRC");

helper_funcs.Plot_RF_Signal(t, real(y_lpf_rrc), "Πραγματικό μέρος σήματος στην έξοδο του LPF με RRC");


%% =========================================================
% 34. Αστερισμός στην έξοδο του LPF
% Κανονικοποίηση μόνο για οπτική σύγκριση
% =========================================================
sym_lpf = y_lpf(delay+1 : sps : end-delay);
sym_lpf_rrc = y_lpf_rrc(delay+1 : sps : end-delay);

sym_lpf_norm = sym_lpf / sqrt(mean(abs(sym_lpf).^2));
sym_lpf_norm_rrc = sym_lpf_rrc / sqrt(mean(abs(sym_lpf_rrc).^2));

helper_funcs.Plot_Constellation(sym_lpf_norm, "Αστερισμός στην έξοδο του LPF χωρίς RRC");

helper_funcs.Plot_Constellation(sym_lpf_norm_rrc, "Αστερισμός στην έξοδο του LPF με RRC");

%% =========================================================
% 35. Υπολογισμός SNR στην έξοδο του LPF
% =========================================================

SNR_lpf_dB = helper_funcs.Calculate_SNR(y_lpf_clean, y_lpf_noise);
fprintf('SNR στην έξοδο του LPF = %.2f dB\n', SNR_lpf_dB);

%% =========================================================
% 36. Υπολογισμός EVM στην έξοδο του LPF
% =========================================================

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_lpf);
fprintf('EVM στην έξοδο του LPF = %.2f %%\n', evm_percent);

%% =========================================================
% 37. 6η βαθμίδα δέκτη: Baseband VGA
% Το VGA αυξάνει τη στάθμη του baseband σήματος πριν από τον ADC.
% Σε αυτό το βήμα βάζουμε μόνο το gain του VGA.
% =========================================================

G_VGA_dB  = 16;                    % gain του VGA από το operating point
G_VGA_lin = 10^(G_VGA_dB/10);

G_BBA_dB = 20;
G_BBA = 10^(G_BBA_dB/10);

% Έξοδος VGA
y_vga = sqrt(G_VGA_lin) * G_BBA * y_lpf_rrc;

%% =========================================================
% 38. Ισχύς στην έξοδο του VGA
% =========================================================
P_vga_out_dBm = 10*log10(mean(abs(y_vga).^2)/1e-3);
fprintf('Ισχύς στην έξοδο του VGA = %.2f dBm\n', P_vga_out_dBm);

P_vga_expected_dBm = P_lpf_out_dBm + G_VGA_dB;
fprintf('Αναμενόμενη ισχύς στην έξοδο του VGA = %.2f dBm\n', P_vga_expected_dBm);

%% =========================================================
% 39. Διάγραμμα σήματος στην έξοδο του VGA
% =========================================================

helper_funcs.Plot_RF_Signal(t, real(y_vga), "Πραγματικό μέρος σήματος στην έξοδο του VGA");

%%
% Second Test
figure;
[pxx, f] = pwelch(y_vga, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Φάσμα Baseband LPF me RRC');
grid on;
xlim([-BW/(2*1e9) - 2  BW/(2*1e9) + 2]);
%% =========================================================
% 40. Αστερισμός στην έξοδο του VGA
% Κανονικοποίηση μόνο για οπτική σύγκριση
% =========================================================
sym_vga = y_vga(delay+1 : sps : end-delay);

sym_vga_norm = sym_vga / sqrt(mean(abs(sym_vga).^2));

helper_funcs.Plot_Constellation(sym_vga_norm, "Αστερισμός στην έξοδο του VGA");

%% =========================================================
% 41. SNR και EVM στην έξοδο του VGA
% Εφόσον εδώ εφαρμόζουμε μόνο gain, το SNR ιδανικά δεν αλλάζει
% =========================================================

y_vga_clean = sqrt(G_VGA_lin) * y_lpf_clean;
y_vga_noise = sqrt(G_VGA_lin) * y_lpf_noise;

SNR_vga_dB = helper_funcs.Calculate_SNR(y_vga_clean, y_vga_noise);
fprintf('SNR στην έξοδο του VGA = %.2f dB\n', SNR_vga_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_vga);
fprintf('EVM στην έξοδο του VGA = %.2f %%\n', evm_percent);

%% =========================================================
% 42. 7η βαθμίδα δέκτη: ADC
% Ο ADC είναι dual, άρα στο complex baseband μοντέλο
% κβαντίζουμε ξεχωριστά το I και το Q κανάλι.
% =========================================================

ADC_bits = 12;                        % ανάλυση ADC
P_FS_ADC_dBm = 3.9;                   % full-scale input power από το datasheet
P_FS_ADC_W = 1e-3 * 10^(P_FS_ADC_dBm/10);

% Για κάθε πραγματικό κανάλι (I ή Q) θεωρούμε full-scale ημιτονοειδές
A_FS_rms  = sqrt(P_FS_ADC_W);         % RMS πλάτος που αντιστοιχεί στο full-scale power
A_FS_peak = sqrt(2) * A_FS_rms;       % peak amplitude για ημίτονο

% Βήμα κβαντοποίησης
Delta = 2*A_FS_peak / (2^ADC_bits - 1);

%% =========================================================
% 43. Διαχωρισμός σε I και Q πριν τον ADC
% =========================================================
I_adc_in = real(y_vga);
Q_adc_in = imag(y_vga);

%% =========================================================
% 44. Clipping στα όρια full-scale
% Αν το σήμα ξεπερνά τα όρια του ADC, κόβεται
% =========================================================
I_clip = min(max(I_adc_in, -A_FS_peak), A_FS_peak);
Q_clip = min(max(Q_adc_in, -A_FS_peak), A_FS_peak);

%% =========================================================
% 45. Κβαντοποίηση 12-bit
% =========================================================
I_adc = Delta * round(I_clip / Delta);
Q_adc = Delta * round(Q_clip / Delta);

% Επανασύνθεση μιγαδικού baseband σήματος
y_adc = I_adc + 1j*Q_adc;

%% =========================================================
% 46. Έλεγχος ισχύος και στάθμης στον ADC
% =========================================================
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
% 47. Διάγραμμα χρόνου: είσοδος και έξοδος ADC (κανάλι I)
% =========================================================

figure;
plot(t*1e9, I_adc_in, 'LineWidth', 1.1);
hold on;
stairs(t*1e9, I_adc, 'LineWidth', 1.1);
grid on;
xlabel('Χρόνος (ns)');
ylabel('Πλάτος');
title('Κανάλι I πριν και μετά τον ADC');
legend('Είσοδος ADC', 'Έξοδος ADC (κβαντισμένη)');
xlim([28 30]);

%% =========================================================
% 48. Αστερισμός στην έξοδο του ADC
% =========================================================
sym_adc = y_adc(delay+1 : sps : end-delay);

sym_adc_norm = sym_adc / sqrt(mean(abs(sym_adc).^2));

helper_funcs.Plot_Constellation(sym_adc_norm, "Αστερισμός στην έξοδο του ADC");

%% =========================================================
% 49. Effective SNR στην έξοδο του ADC
% Το σφάλμα τώρα περιλαμβάνει:
% - αναλογικό θόρυβο προηγούμενων βαθμίδων
% - κβαντοποίηση
% - πιθανό clipping
% =========================================================

SNR_adc_eff_dB = helper_funcs.Calculate_SNR(y_vga_clean, y_adc - y_vga_clean);
fprintf('Effective SNR στην έξοδο του ADC = %.2f dB\n', SNR_adc_eff_dB);

%% =========================================================
% 50. EVM στην έξοδο του ADC
% =========================================================

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_adc);
fprintf('EVM στην έξοδο του ADC = %.2f %%\n', evm_percent);