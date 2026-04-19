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
noise_in = sqrt(Nin * 50 /2) * ...
    (randn(size(rx_in)) + 1j*randn(size(rx_in)));

rx_in_noisy = rx_in + noise_in;

P_rx_check_dBm = 10*log10(mean(abs(rx_in).^2)/1e-3/50);
fprintf('Ισχύς στην είσοδο του δέκτη = %.2f dBm\n', P_rx_check_dBm);

%% =========================================================
% 15. 2η βαθμίδα δέκτη: LNA

G_LNA_dB  = 25;               
G_LNA_lin = 10^(G_LNA_dB/10);
NF_LNA_dB = 5;
NF_LNA = helper_funcs.dB_To_Linear(NF_LNA_dB, 'dB');

% Έξοδος LNA στο complex baseband:
% το πλάτος του σήματος πολλαπλασιάζεται με sqrt(G)
y_lna_clean = sqrt(G_LNA_lin) * rx_in;
noise_after_lna_gain = sqrt(G_LNA_lin) * noise_in;
y_lna = sqrt(G_LNA_lin) * rx_in_noisy;
SNR_out_LNA = SNRin_dB - NF_LNA_dB;
[y_lna_noisy, LNA_NF_noise] = helper_funcs.Calcualte_noisy_signal(y_lna, SNRin_dB, SNR_out_LNA);
noise_after_LNA = noise_after_lna_gain + LNA_NF_noise;

% Real Part of the signal after the LNA with noise
helper_funcs.Plot_RF_Signal(t, real(y_lna_noisy), "Real part of signal after LNA");

% Constellation after the LNA with noise
sym_lna = y_lna_noisy(delay+1 : sps : end-delay);

sym_lna_noisy_norm = sym_lna / sqrt(mean(abs(sym_lna).^2));

helper_funcs.Plot_Constellation(sym_lna_noisy_norm, "Αστερισμός στην έξοδο του LNA");

% Power calculations after LNA
% We expect P_LNA,out = P_BPF,out + 25 dB
P_lna_out_dBm = 10*log10(mean(abs(y_lna_noisy).^2)/1e-3/50);
P_lna_expected_dBm = PRX_dBm + G_LNA_dB;

fprintf('Expected Power after LNA: %.2f dBm\n', P_lna_expected_dBm);
fprintf('Power after LNA: %.2f dBm\n', P_lna_out_dBm);

% Now we calculate the SNR and EVM after LNA

SNR_lna_meas_dB = helper_funcs.Calculate_SNR(y_lna_clean, noise_after_LNA);
fprintf('SNR στην έξοδο του LNA = %.2f dB\n', SNR_lna_meas_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_lna);
fprintf('EVM στην έξοδο του LNA = %.2f %%\n', evm_percent);

%% ========================================================= 
% 1η βαθμίδα δέκτη: RF BPF ως ισοδύναμη baseband επίδραση
G_BPF_dB = -1.5;                          % insertion loss 1.5 dB
G_BPF_lin = 10^(G_BPF_dB/10);
NF_BPF_dB = - G_BPF_dB; % Not mentioned in the datasheet
NF_BPF = helper_funcs.dB_To_Linear(NF_BPF_dB, 'dB');

y_bpf_clean = sqrt(G_BPF_lin) * y_lna_clean;
noise_after_bpf_gain = sqrt(G_BPF_lin) * noise_after_LNA; 
y_bpf = sqrt(G_BPF_lin) * y_lna_noisy;
F_BPF = NF_LNA + (NF_BPF - 1) / G_LNA_lin;
SNR_out_BPF_dB = SNR_out_LNA - 10 * log10(F_BPF);
[y_bpf_noisy, bpf_NF_noise] = helper_funcs.Calcualte_noisy_signal(y_bpf, SNR_out_LNA, SNR_out_BPF_dB);
noise_after_bpf = noise_after_bpf_gain + bpf_NF_noise;

P_bpf_out_dBm = 10*log10(mean(abs(y_bpf_noisy).^2)/1e-3/50);
P_bpf_expected_dBm = P_lna_expected_dBm + G_BPF_dB;
fprintf('Expected Power after BPF: %.2f dBm\n', P_bpf_expected_dBm);
fprintf('Power after BPF: %.2f dBm\n', P_bpf_out_dBm);

% Real Part of the signal after the BPF with noise
helper_funcs.Plot_RF_Signal(t, real(y_bpf_noisy), "Πραγματικό μέρος σήματος στην έξοδο του RF BPF");

% Constellation after the BPF with noise
sym_bpf = y_bpf_noisy(delay+1 : sps : end-delay);

% We normalize so that we can properly compare
sym_bpf_norm = sym_bpf / sqrt(mean(abs(sym_bpf).^2));

helper_funcs.Plot_Constellation(sym_bpf_norm, "Αστερισμός στην έξοδο του RF BPF");

SNR_out_BPF_calc_dB = helper_funcs.Calculate_SNR(y_bpf, noise_after_bpf);
fprintf('SNR στην έξοδο του BPF = %.2f dB\n', SNR_out_BPF_calc_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_bpf);
fprintf('EVM στην έξοδο του BPF = %.2f %%\n', evm_percent);

%% =========================================================
% 3η βαθμίδα δέκτη: RF Splitter


G_PS_dB  = -3.5;                    % συνολική απώλεια ανά κλάδο
G_PS_lin = 10^(G_PS_dB/10);
NF_PS_dB = - G_PS_dB; % Not mentioned in the datasheet
NF_PS = helper_funcs.dB_To_Linear(NF_PS_dB, 'dB');

% Έξοδοι splitter στους δύο κλάδους
y_I_split_clean = sqrt(G_PS_lin) * y_bpf_clean;
noise_after_split_gain = sqrt(G_PS_lin) * noise_after_bpf;
y_I_split = sqrt(G_PS_lin) * y_bpf_noisy;
y_Q_split = sqrt(G_PS_lin) * y_bpf_noisy;
F_I_Split = F_BPF + (NF_PS - 1) / (G_BPF_lin * G_LNA_lin);
F_Q_Split = F_BPF + (NF_PS - 1) / (G_BPF_lin * G_LNA_lin);
SNR_out_I_split = SNR_out_LNA - 10 * log10(F_I_Split);
SNR_out_Q_split = SNR_out_LNA - 10 * log10(F_Q_Split);
[y_I_split_noisy, split_NF_I_noise] = helper_funcs.Calcualte_noisy_signal(y_I_split, SNR_out_LNA, SNR_out_I_split);
[y_Q_split_noisy, split_NF_Q_noise] = helper_funcs.Calcualte_noisy_signal(y_Q_split, SNR_out_LNA, SNR_out_Q_split);
noise_after_split = noise_after_split_gain + split_NF_I_noise;

helper_funcs.Plot_RF_Signal(t, real(y_I_split_noisy), "Real part of signal after Power Splitter (Branch I)");

sym_I_split = y_I_split_noisy(delay+1 : sps : end-delay);

sym_I_split_norm = sym_I_split / sqrt(mean(abs(sym_I_split).^2));

helper_funcs.Plot_Constellation(sym_I_split_norm, "Constellation after RF Splitter (branch I)");

% Υπολογισμός ισχύος ανά κλάδο
P_I_split_dBm = 10*log10(mean(abs(y_I_split_noisy).^2) /1e-3 / 50);
P_Q_split_dBm = 10*log10(mean(abs(y_Q_split_noisy).^2) /1e-3 / 50);
P_split_expected_dBm = P_bpf_expected_dBm + G_PS_dB;

fprintf('Expected Power per branch after Power Splitter: %.2f dBm\n', P_split_expected_dBm);
fprintf('Power at branch I after Power Splitter: %.2f dBm\n', P_I_split_dBm);
fprintf('Power at branch Q after Power Splitter: %.2f dBm\n', P_Q_split_dBm);

% Now we calculate the SNR and EVM after Power Spliiter

SNR_PS_dB = helper_funcs.Calculate_SNR(y_I_split_clean, noise_after_split);
fprintf('SNR στην έξοδο του RF Splitter = %.2f dB\n', SNR_PS_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_I_split);
fprintf('EVM στην έξοδο του RF Splitter = %.2f %%\n', evm_percent);

%% =========================================================
% 4η βαθμίδα δέκτη: Mixer

G_MixerI_dB  = -9;                    % conversion loss από το datasheet
G_MixerI_lin = 10^(G_MixerI_dB/10);
NF_MixerI_dB = - G_MixerI_dB;
NF_MixerI = helper_funcs.dB_To_Linear(NF_MixerI_dB, 'dB');

% We use just one branch of the split
y_mix_clean = sqrt(G_MixerI_lin) * y_I_split_clean;
y_mix = sqrt(G_MixerI_lin) * y_I_split_noisy;
noise_after_mix_gain = sqrt(G_MixerI_lin) * noise_after_split;
F_Mixer = F_I_Split + (NF_MixerI - 1) / (G_BPF_lin * G_LNA_lin * G_PS_lin);
SNR_out_mixer_dB = SNR_out_I_split - 10 * log10(F_Mixer);
[y_mix_noisy, mixer_NF_noise] = helper_funcs.Calcualte_noisy_signal(y_mix, SNR_out_I_split, SNR_out_mixer_dB);
noise_after_mix = noise_after_mix_gain + mixer_NF_noise;

helper_funcs.Plot_RF_Signal(t, real(y_mix_noisy), "Real Part of Signal after Mixer");

sym_mix = y_mix(delay+1 : sps : end-delay);

sym_mix_norm = sym_mix / sqrt(mean(abs(sym_mix).^2));

helper_funcs.Plot_Constellation(sym_mix_norm, "Constellation after Mixer");

% Υπολογισμός ισχύος στην έξοδο του mixer
P_mix_out_dBm = 10*log10(mean(abs(y_mix_noisy).^2)/1e-3/50);
P_mix_expected_dBm = P_split_expected_dBm + G_MixerI_dB;

fprintf('Αναμενόμενη ισχύς στην έξοδο του Mixer = %.2f dBm\n', P_mix_expected_dBm);
fprintf('Ισχύς στην έξοδο του Mixer = %.2f dBm\n', P_mix_out_dBm);

SNR_mix_dB = helper_funcs.Calculate_SNR(y_mix_clean, noise_after_mix);
fprintf('SNR στην έξοδο του Mixer = %.2f dB\n', SNR_mix_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_mix);
fprintf('EVM στην έξοδο του Mixer = %.2f %%\n', evm_percent);

%% =========================================================
% 5η βαθμίδα δέκτη: LPF and RRC

G_LPF_dB  = -3.1;                   % insertion loss του LPF
G_LPF_lin = 10^(G_LPF_dB/10);
NF_LPF_dB = - G_LPF_dB; % Not mentioned in the datasheet
NF_LPF = helper_funcs.dB_To_Linear(NF_LPF_dB, 'dB');


fc_lpf = 3e9;                       % cutoff συχνότητα ~ 3 GHz
Wn = fc_lpf / (Fs/2);               % κανονικοποιημένη συχνότητα για fir1

lpfOrder = 80;                      % τάξη FIR φίλτρου
b_lpf = fir1(lpfOrder, Wn, 'low');  % συντελεστές LPF
b_lpf = b_lpf / sum(b_lpf);

% LPF
y_lpf_clean = sqrt(G_LPF_lin) * conv(y_mix_clean, b_lpf, 'same');
y_lpf_noisy = sqrt(G_LPF_lin) * conv(y_mix_noisy, b_lpf, 'same');
noise_after_lpf = sqrt(G_LPF_lin) * conv(noise_after_mix, b_lpf, 'same');

% RRC
% The RRC does not add noise through the Noise Figure as it is normally implemented at a DSP level
y_rrc_clean = conv(y_lpf_clean, rrc, 'same');
y_rrc_noisy = conv(y_lpf_noisy, rrc, 'same'); 
noise_after_rrc = conv(noise_after_lpf, b_lpf, 'same');

% Power Normalization to ensure the RRC does not affect the power

% Clean Signal
power_before_clean = mean(abs(y_lpf_clean).^2);
power_after_clean = mean(abs(y_rrc_clean).^2);
y_rrc_clean = y_rrc_clean * sqrt(power_before_clean / power_after_clean);

% Noisy Signal
power_before_noisy = mean(abs(y_lpf_noisy).^2);
power_after_noisy = mean(abs(y_rrc_noisy).^2);
y_rrc_noisy = y_rrc_noisy * sqrt(power_before_noisy / power_after_noisy);

% Noise
noise_power_before_rrc = mean(abs(noise_after_lpf).^2);
noise_power_after_rrc = mean(abs(noise_after_rrc).^2);
noise_after_rrc = noise_after_rrc * sqrt(noise_power_before_rrc / noise_power_after_rrc);

helper_funcs.Plot_RF_Signal(t, real(y_rrc_noisy), "LPF Signal with RRC");

sym_rrc = y_rrc_noisy(delay+1 : sps : end-delay);

sym_lpf_norm = sym_rrc / sqrt(mean(abs(sym_rrc).^2));

helper_funcs.Plot_Constellation(sym_lpf_norm, "Constellation after LPF and RRC");

P_rrc_out_dBm = 10*log10(mean(abs(y_rrc_noisy).^2)/1e-3/50);
P_rrc_expected_dBm = P_mix_expected_dBm + G_LPF_dB;

fprintf('Expected Power after LPF and RRC: %.2f dBm\n', P_rrc_expected_dBm);
fprintf('Power after LPF and RRC: %.2f dBm\n', P_rrc_out_dBm);

SNR_out_rrc_dB = helper_funcs.Calculate_SNR(y_rrc_clean, noise_after_rrc);
fprintf('SNR στην έξοδο του LPF και RRC = %.2f dB\n', SNR_out_rrc_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_rrc);
fprintf('EVM στην έξοδο του LPF = %.2f %%\n', evm_percent);


% Effect of LPF RRC at the frequency domain
figure;
[pxx, f] = pwelch(y_mix_noisy, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Frequency (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Baseband spectre after mixer');
grid on;
xlim([-BW/(2*1e9) - 2  BW/(2*1e9) + 2]);

figure;
[pxx, f] = pwelch(y_lpf_noisy, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Baseband spectre after LPF without RRC');
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
[pxx, f] = pwelch(y_rrc_noisy, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Φάσμα Baseband LPF me RRC');
grid on;
xlim([-BW/(2*1e9) - 2  BW/(2*1e9) + 2]);

%%
% 6η βαθμίδα δέκτη: Cascade Baseband Amplifiers

G_BBA_dB = 20;
G_BBA_lin = helper_funcs.dB_To_Linear(G_BBA_dB, 'dB');
NF_BBA_dB = 12.1;
NF_BBA = helper_funcs.dB_To_Linear(NF_BBA_dB, 'dB');

% Stage 1
y_BBA1_clean = sqrt(G_BBA_lin) * y_rrc_clean;
y_BBA1 = sqrt(G_BBA_lin) * y_rrc_noisy;
noise_after_BBA1_gain = sqrt(G_BBA_lin) * noise_after_rrc;
F_LPF = 10 ^ ((SNR_out_mixer_dB - SNR_out_rrc_dB) / 10);
F_BBA1 = F_LPF + (NF_BBA - 1) / (G_BPF_lin * G_LNA_lin * G_PS_lin * G_MixerI_lin * G_LPF_lin);
SNR_out_BBA1_dB = SNR_out_rrc_dB - 10 * log10(F_BBA1);
[y_BBA1_noisy, BBA1_NF_noise] = helper_funcs.Calcualte_noisy_signal(y_BBA1, SNR_out_rrc_dB, SNR_out_BBA1_dB);
noise_after_BBA1 = noise_after_BBA1_gain + BBA1_NF_noise;

% Stage 2
y_BBA2_clean = sqrt(G_BBA_lin) * y_BBA1_clean;
y_BBA2 = sqrt(G_BBA_lin) * y_BBA1_noisy;
noise_after_BBA2_gain = sqrt(G_BBA_lin) * noise_after_BBA1;
F_BBA2 = F_BBA1 + (NF_BBA - 1) / (G_BPF_lin * G_LNA_lin * G_PS_lin * G_MixerI_lin * G_LPF_lin * G_BBA_lin);
SNR_out_BBA2_dB = SNR_out_BBA1_dB - 10 * log10(F_BBA2);
[y_BBA2_noisy, BBA2_NF_noise] = helper_funcs.Calcualte_noisy_signal(y_BBA2, SNR_out_rrc_dB, SNR_out_BBA2_dB);
noise_after_BBA2 = noise_after_BBA2_gain + BBA2_NF_noise;

helper_funcs.Plot_RF_Signal(t, real(y_BBA2_noisy), "Real part of the signal after the BBAs");

sym_BBA = y_BBA2_noisy(delay+1 : sps : end-delay);

sym_BBA_norm = sym_BBA / sqrt(mean(abs(sym_BBA).^2));

helper_funcs.Plot_Constellation(sym_BBA_norm, "Constellation after LPF and RRC");

P_BBA_out_dBm = 10*log10(mean(abs(y_BBA2_noisy).^2)/1e-3/50);
P_BBA_expected_dBm = P_rrc_expected_dBm + 2 * G_BBA_dB;

fprintf('Expected Power after the BBAs: %.2f dBm\n', P_BBA_expected_dBm);
fprintf('Power after the BBAs: %.2f dBm\n', P_BBA_out_dBm);

SNR_BBA_dB = helper_funcs.Calculate_SNR(y_BBA2_clean, noise_after_BBA2);
fprintf('SNR στην έξοδο του LPF και RRC = %.2f dB\n', SNR_BBA_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_BBA);
fprintf('EVM στην έξοδο του LPF = %.2f %%\n', evm_percent);

%%
% Second LPF

y_lpf2_clean = sqrt(G_LPF_lin) * conv(y_BBA2_clean, b_lpf, 'same');
y_lpf2_noisy = sqrt(G_LPF_lin) * conv(y_BBA2_noisy, b_lpf, 'same');
noise_after_lpf2 = sqrt(G_LPF_lin) * conv(noise_after_BBA2, b_lpf, 'same');

sym_lpf2 = y_lpf2_noisy(delay+1 : sps : end-delay);

sym_lpf2_norm = sym_rrc / sqrt(mean(abs(sym_lpf2).^2));


P_lpf2_out_dBm = 10*log10(mean(abs(y_lpf2_noisy).^2)/1e-3/50);
P_lpf2_expected_dBm = P_BBA_expected_dBm + G_LPF_dB;

fprintf('Expected Power after LPF2: %.2f dBm\n', P_lpf2_expected_dBm);
fprintf('Power after LPF2: %.2f dBm\n', P_lpf2_out_dBm);

SNR_out_lpf2_dB = helper_funcs.Calculate_SNR(y_lpf2_clean, noise_after_lpf2);
fprintf('SNR στην έξοδο του LPF 2: %.2f dB\n', SNR_out_lpf2_dB);

evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_lpf2);
fprintf('EVM στην έξοδο του LPF = %.2f %%\n', evm_percent);

%% Testing a third BBA
y_BBA3_clean = sqrt(G_BBA_lin) * y_lpf2_clean;
y_BBA3 = sqrt(G_BBA_lin) * y_lpf2_noisy;
noise_after_BBA3_gain = sqrt(G_BBA_lin) * noise_after_lpf2;
F_LPF2 = 10 ^ ((SNR_out_BBA2_dB - SNR_out_lpf2_dB) / 10);
F_BBA3 = F_LPF2 + (NF_BBA - 1) / (G_BPF_lin * G_LNA_lin * G_PS_lin * G_MixerI_lin * G_LPF_lin^2 * G_BBA_lin ^2);
SNR_out_BBA3_dB = SNR_out_lpf2_dB - 10 * log10(F_BBA3);
[y_BBA3_noisy, BBA3_NF_noise] = helper_funcs.Calcualte_noisy_signal(y_BBA3, SNR_out_lpf2_dB, SNR_out_BBA3_dB);
noise_after_BBA3 = noise_after_BBA3_gain + BBA3_NF_noise;



%% =========================================================
% 37. 6η βαθμίδα δέκτη: Baseband VGA

G_VGA_dB  = 0;                    % gain του VGA από το operating point
G_VGA_lin = 10^(G_VGA_dB/10);
NF_VGA_dB = 14.7;
NF_VGA_lin = helper_funcs.dB_To_Linear(NF_VGA_dB, 'dB');

y_vga_clean = sqrt(G_VGA_lin) * y_BBA3_clean;
y_vga = sqrt(G_VGA_lin) * y_BBA3_noisy;
noise_after_vga_gain = sqrt(G_VGA_lin) * noise_after_BBA3;
%F_LPF2 = 10 ^ ((SNR_out_BBA2_dB - SNR_out_lpf2_dB) / 10);
%F_vga = F_LPF2 + (NF_VGA_lin - 1) / (G_BPF_lin * G_LNA_lin * G_PS_lin * G_MixerI_lin * G_LPF_lin^2 * G_BBA_lin ^2);
F_vga = F_BBA3 + (NF_VGA_lin - 1) / (G_BPF_lin * G_LNA_lin * G_PS_lin * G_MixerI_lin * G_LPF_lin^2 * G_BBA_lin ^3);
SNR_out_vga_dB = SNR_out_BBA3_dB - 10 * log10(F_vga);
[y_vga_noisy, vga_NF_noise] = helper_funcs.Calcualte_noisy_signal(y_vga, SNR_out_lpf2_dB, SNR_out_vga_dB);
noise_after_vga = noise_after_vga_gain + vga_NF_noise;

helper_funcs.Plot_RF_Signal(t, real(y_vga_noisy), "Real part of signal after vga");

sym_vga = y_vga_noisy(delay+1 : sps : end-delay);

sym_vga_norm = sym_vga / sqrt(mean(abs(sym_vga).^2));

helper_funcs.Plot_Constellation(sym_vga_norm, "Αστερισμός στην έξοδο του VGA");


P_vga_out_dBm = 10*log10(mean(abs(y_vga_noisy).^2)/1e-3/50);
P_vga_expected_dBm = P_lpf2_expected_dBm + G_VGA_dB + G_BBA_dB;

fprintf('Expected Power after VGA: %.2f dBm\n', P_vga_expected_dBm);
fprintf('Power after VGA: %.2f dBm\n', P_vga_out_dBm);

SNR_out_vga_dB = helper_funcs.Calculate_SNR(y_vga_clean, noise_after_vga);
fprintf('SNR στην έξοδο του VGA: %.2f dB\n', SNR_out_vga_dB);


evm_percent = helper_funcs.Calculate_EVM(sym_samples, sym_vga);
fprintf('EVM στην έξοδο του VGA: %.2f %%\n', evm_percent);

%% =========================================================
% 42. 7η βαθμίδα δέκτη: ADC
% Ο ADC είναι dual, άρα στο complex baseband μοντέλο
% κβαντίζουμε ξεχωριστά το I και το Q κανάλι.
% =========================================================

ADC_bits = 12;                        % ανάλυση ADC
P_FS_ADC_dBm = 3.9;                   % full-scale input power από το datasheet
P_FS_ADC_W = 1e-3 * 10^(P_FS_ADC_dBm/10);

% Για κάθε πραγματικό κανάλι (I ή Q) θεωρούμε full-scale ημιτονοειδές
A_FS_rms  = sqrt(P_FS_ADC_W * 50);         % RMS πλάτος που αντιστοιχεί στο full-scale power
A_FS_peak = sqrt(2) * A_FS_rms;       % peak amplitude για ημίτονο

% Βήμα κβαντοποίησης
Delta = 2*A_FS_peak / (2^ADC_bits - 1);

%% =========================================================
% 43. Διαχωρισμός σε I και Q πριν τον ADC
% =========================================================
I_adc_in = real(y_vga_noisy);
Q_adc_in = imag(y_vga_noisy);

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