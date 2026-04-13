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
filtlen = 20; % Filter cut-off at 20 symbols
Rs = BW/(1 + rolloff);
Fs = sps * Rs; %συχνοτητα δειγματοληψιας, πρεπει τουλαχιστον 2 επι το symbol rate
N_syms = 10000;
bits = randi([0,1], 1 , N_syms * k);
idx = bi2de(reshape(bits, k , []).' , 'left-msb');
% Rotate PSK by pi/4 so the values are not on the axis
syms = pskmod(idx, M, pi/4);
% filtlen*sps + 1 = 81 so 81 coefficients, odd because rc is
% symmetrical around 0
rrc = rcosdesign( rolloff, filtlen, sps, 'sqrt');
% umpsample so each symbol becomes sps = 8 samples
% rrc in each symbol
tx_shaped = upfirdn(syms, rrc, sps);
delay = (filtlen * sps) / 2;
% Keep only after the transient delay introduced by the filter
tx_shaped_valid = tx_shaped(delay + 1: end - delay);
% Take the I and Q of the complex values
tx_I = real(tx_shaped_valid');
tx_Q = imag(tx_shaped_valid');
N = length(tx_shaped_valid); % total number of samples
t = (0:N-1)/Fs; % Time between consecutive samples
%% PLOT 1: Χρόνος Baseband
figure;
subplot(2,1,1);
plot(t*1e9, tx_I, 'b', 'LineWidth', 1.2);
xlabel('Χρόνος (ns)'); ylabel('Πλάτος');
title('QPSK Baseband — I'); grid on;
xlim([0, 100]);
subplot(2,1,2);
plot(t*1e9, tx_Q, 'r', 'LineWidth', 1.2);
xlabel('Χρόνος (ns)'); ylabel('Πλάτος');
title('QPSK Baseband — Q'); grid on;
xlim([0, 100]);
% Keep just the first 100 nanoseconds to see the value more clearly
%% PLOT 2: Φάσμα Baseband
figure;
[pxx, f] = pwelch(tx_shaped_valid, [], [], [], Fs/1e9, 'centered');
plot(f, 10*log10(pxx));
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
title('Φάσμα Baseband QPSK');
grid on;
xlim([-BW/(2*1e9) - 2  BW/(2*1e9) + 2]);
% Add +/- 2GHz so we can see the RRC filter
%% PLOT 3: Αστερισμός
figure;
plot(real(syms), imag(syms), 'b.', 'MarkerSize', 15);
xlabel('I'); ylabel('Q');
title('Αστερισμός QPSK');
grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);
%% PLOT 4: Φάσμα RF
tx_passband = real(tx_shaped_valid .* exp(1j * 2 * pi * fc * t'));
figure;
[pxx_pb, f_pb] = pwelch(tx_passband, [], [], [], Fs/1e9);
plot(f_pb, 10 * log(pxx_pb));
title('Φάσμα Passband QPSK');
xlabel('Συχνότητα (GHz)');
ylabel('Power/frequency (dB/Hz)');
grid on;
xlim([fc/1e9 - 5, fc/1e9 + 5]);
%Same here as the baseband
%% PLOT 5: RF - Χρόνος
figure;
% Plot the signal, converting the time vector to nanoseconds (ns) for readability
plot(t' * 1e9, tx_passband);
title('Passband Signal in Time Domain (RF)');
xlabel('Time (ns)');
ylabel('Amplitude');
grid on;

% Zoom in to see the actual carrier waves
% Here we are looking at just the first 2 nanoseconds of the transmission.
xlim([0, 2]);
%% PLOT 6: Αστερισμός μετά από RRC
figure;
plot(real(tx_shaped_valid(1000: 2000 + 1)), imag(tx_shaped_valid(1000: 2000 + 1)), 'b.', 'MarkerSize', 15);
xlabel('I'); ylabel('Q');
title('Αστερισμός QPSK μετά από το RRC φίλτρο');
grid on; axis equal;
xlim([-1 1]); ylim([-1 1]);
%% PLOT 8: Eye Diagram
% figure;
% eyediagram(tx_I, 2*sps);
% title('Eye Diagram QPSK');
% %% PLOT 8: Απόκριση RRC φίλτρου
% figure;
% [H, f] = freqz(rrc, 1, 1024, Fs/1e9);
% subplot(2,1,1);
% plot(f, 20*log10(abs(H)), 'b', 'LineWidth', 1.5);
% xlabel('Συχνότητα (GHz)');
% ylabel('Πλάτος (dB)');
% title('RRC Φίλτρο — Απόκριση Συχνότητας');
% grid on;
% xlim([0 BW/1e9]);
% subplot(2,1,2);
% plot(f, unwrap(angle(H))*180/pi, 'r', 'LineWidth', 1.5);
% xlabel('Συχνότητα (GHz)');
% ylabel('Φάση (°)');
% title('RRC Φίλτρο — Απόκριση Φάσης');
% grid on;
% xlim([0 BW/1e9]);
% %% PLOT: RRC παλμός στο χρόνο
% figure;
% t_rrc = (-filtlen*sps/2 : filtlen*sps/2) / Fs * 1e9;
% plot(t_rrc, rrc, 'b', 'LineWidth', 1.5);
% xlabel('Χρόνος (ns)');
% ylabel('Πλάτος');
% title('Root Raised Cosine Παλμός');
% grid on;
% xline(0, 'r--', 'LineWidth', 1.5);  % κέντρο συμμετρίας


