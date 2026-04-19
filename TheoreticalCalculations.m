clear;
clc;
close all

%%
% RF Bandpass Filter (BPF) Parameters
G_BPF_dB = -1.5;
G_BPF = helper_funcs.dB_To_Linear(G_BPF_dB, 'dB');
NF_BPF_dB = - G_BPF_dB; % Not mentioned in the datasheet
NF_BPF = helper_funcs.dB_To_Linear(NF_BPF_dB, 'dB');

%%
% LNA Parameters
G_LNA_dB = 25;
G_LNA = helper_funcs.dB_To_Linear(G_LNA_dB, 'dB');
NF_LNA_dB = 5;
NF_LNA = helper_funcs.dB_To_Linear(NF_LNA_dB, 'dB');
P1dB_LNA_dBm = -3;
OIP3_LNA_dBm = P1dB_LNA_dBm + 10;
OIP3_LNA = helper_funcs.dB_To_Linear(OIP3_LNA_dBm, 'dBm');
IIP3_LNA_dBm = P1dB_LNA_dBm + 10 - G_LNA_dB;

%%
% RF Power Splitter (PS) Parameters
G_PS_dB = -3.5;
G_PS = helper_funcs.dB_To_Linear(G_PS_dB, 'dB');
NF_PS_dB = - G_PS_dB; % Not mentioned in the datasheet
NF_PS = helper_funcs.dB_To_Linear(NF_PS_dB, 'dB');
Isolation_PS_dB = 20;

%%
% Mixer I
G_MixerI_dB = -9;
G_MixerI = helper_funcs.dB_To_Linear(G_MixerI_dB, 'dB');
NF_MixerI_dB = - G_MixerI_dB;
NF_MixerI = helper_funcs.dB_To_Linear(NF_MixerI_dB, 'dB');
P1dB_MixerI_dBm = -3;
Isolation_MixerI_dB = 30;
LO_Power_MixerI_dBm = 13;
% Since we have a passive mixer
IIP3_MixerI_dBm = LO_Power_MixerI_dBm + 10;
OIP3_MixerI_dBm = IIP3_MixerI_dBm + G_MixerI_dB;
OIP3_MixerI = helper_funcs.dB_To_Linear(OIP3_MixerI_dBm, 'dBm');

%%
% Mixer Q
G_MixerQ_dB = G_MixerI_dB;
G_MixerQ = G_MixerI;
NF_MixerQ_dB = NF_MixerI_dB;
NF_MixerQ = NF_MixerI;
P1dB_MixerQ_dBm = P1dB_MixerI_dBm;
Isolation_MixerQ_dB = Isolation_MixerI_dB;
LO_Power_MixerQ_dBm = 13;
% Since we have a passive mixer
IIP3_MixerQ_dB = IIP3_MixerI_dBm;
OIP3_MixerQ_dB = OIP3_MixerI_dBm;
OIP3_MixerQ = OIP3_MixerI;

%%
% I-channel LPF (ILPF)
G_ILPF_dB = -3.1;
G_ILPF = helper_funcs.dB_To_Linear(G_ILPF_dB, 'dB');
NF_ILPF_dB = - G_ILPF_dB; % Not mentioned in the datasheet
NF_ILPF = helper_funcs.dB_To_Linear(NF_ILPF_dB, 'dB');

%%
% Q-channel LPF (QLPF)
G_QLPF_dB = G_ILPF_dB;
G_QLPF = G_ILPF;
NF_QLPF_dB = NF_ILPF_dB;
NF_QLPF = NF_ILPF;

%%
% Baseband Amplifier (BBA)
G_BBA_dB = 20;
G_BBA = helper_funcs.dB_To_Linear(G_BBA_dB, 'dB');
NF_BBA_dB = 12.1;
NF_BBA = helper_funcs.dB_To_Linear(NF_BBA_dB, 'dB');
P1dB_BBA_dBm = 14.3;
OIP3_BBA_dBm = 25;
OIP3_BBA = helper_funcs.dB_To_Linear(OIP3_BBA_dBm, 'dBm');
IIP3_BBA_dBm = OIP3_BBA_dBm - G_BBA_dB;
%%
% I-channel Baseband VGA (IVGA)
G_IVGA_dB = 0;
G_IVGA = helper_funcs.dB_To_Linear(G_IVGA_dB, 'dB');
NF_IVGA_dB = 14.7;
NF_IVGA = helper_funcs.dB_To_Linear(NF_IVGA_dB, 'dB');
P1dB_IVGA_dBm = -0.9;
OIP3_IVGA_dBm = 17;
OIP3_IVGA = helper_funcs.dB_To_Linear(OIP3_IVGA_dBm, 'dBm');
IIP3_IVGA_dBm = OIP3_IVGA_dBm - G_IVGA_dB;
%%
% Q-channel Baseband VGA (QVGA)
G_QVGA_dB = G_IVGA_dB;
G_QVGA = G_IVGA;
NF_QVGA_dB = NF_IVGA_dB;
NF_QVGA = NF_IVGA;
OIP3_QVGA = OIP3_IVGA;
P1dB_QVGA_dBm = P1dB_IVGA_dBm;
OIP3_QVGA_dBm = OIP3_IVGA_dBm;
IIP3_QVGA_dBm = IIP3_IVGA_dBm;


%%
% ADC
ADC_Resolution = 12;
NF_ADC_dB = 25.3;
NF_ADC = helper_funcs.dB_To_Linear(NF_ADC_dB, 'dB');

%%
% Initial Calculations
SNRin_dB = 20; % SNR at the receiver
k = 1.38e-23;  % Boltzman constant
T_A = 150;     % In Kelvin
B = 6e+9;      % Bandwidth

Nin = k * T_A * B; % The formula for the noise

% We know that SNR = 10log(S/N) so S = (10 ^ (SNR / 10)) * N
Sin = (10 ^ (SNRin_dB / 10)) * Nin;
Sin_dBm = 10 * log10(Sin / 1e-3);

%%
% Total Gain Calculation
G_Total = G_BPF_dB + G_LNA_dB + G_PS_dB + G_MixerI_dB + 2 * G_ILPF_dB + 3 * G_BBA_dB + G_IVGA_dB;

Voltage_final = 705e-6 * 10 ^ (G_Total / 20);

fprintf("Total gain: %f dB\nFinal Vpp: %f V\n", G_Total, Voltage_final);

%%
% Total Noise Figure Calculation using 
% For passive elements it's F = 1 + (L - 1) * T / T0, where L = 1 / G.
% To calculate the NF in this cascade connection, we use Friis' Formula
NF_total = NF_LNA +...
    (NF_BPF - 1) / G_LNA +...
    (NF_PS - 1) / (G_BPF * G_LNA) +...
    (NF_MixerI - 1) / (G_BPF * G_LNA * G_PS) +...
    (NF_ILPF - 1) / (G_BPF * G_LNA * G_PS * G_MixerI) +...
    (NF_BBA - 1) / (G_BPF * G_LNA * G_PS * G_MixerI * G_ILPF) +...
    (NF_BBA - 1) / (G_BPF * G_LNA * G_PS * G_MixerI * G_ILPF * G_BBA) +...
    (NF_ILPF - 1) / (G_BPF * G_LNA * G_PS * G_MixerI * G_ILPF * G_BBA * G_BBA) +...
    (NF_BBA - 1) / (G_BPF * G_LNA * G_PS * G_MixerI * G_ILPF * G_BBA * G_BBA * G_ILPF) +...
    (NF_IVGA - 1) / (G_BPF * G_LNA * G_PS * G_MixerI * G_ILPF * G_BBA * G_BBA * G_ILPF * G_BBA);

NF_total_dB = helper_funcs.Linear_To_dB(NF_total, 'dB');

fprintf("Total Noise Factor: %f\nTotal Noise Figure: %f dB\n", NF_total, NF_total_dB);

%%
% Total IP3 calculation
OIP3_total_linear_reciprocal =  1 / OIP3_IVGA + 1 / (G_IVGA * OIP3_BBA) +...
    1 / (G_IVGA * G_BBA * OIP3_BBA) + 1 / (G_IVGA * G_BBA * G_BBA * OIP3_MixerI) +...
    1 / (G_IVGA * G_BBA * G_BBA * G_MixerI * OIP3_LNA);
OIP3_total = 1 / OIP3_total_linear_reciprocal;
OIP3_total_dBm = helper_funcs.Linear_To_dB(OIP3_total, 'dBm');
fprintf("Total OIP3: %f dBm\n", OIP3_total_dBm);
%%
% Dynamic Range calculation
MDS = -174 + 10 * log10(B / 2) + NF_total_dB;
P1dB_Input_VGA = P1dB_IVGA_dBm - 2 * helper_funcs.Linear_To_dB(G_Total, 'dB') + G_IVGA_dB;
P1dB_Input_BBA2 = P1dB_BBA_dBm - 2 * helper_funcs.Linear_To_dB(G_Total, 'dB') + G_IVGA_dB + G_BBA_dB;
P1dB_Input_BB1 = P1dB_BBA_dBm - 2 * helper_funcs.Linear_To_dB(G_Total, 'dB') + G_IVGA_dB + 2 * G_BBA_dB;
P1dB_Input_MixerI = P1dB_MixerI_dBm - G_BPF_dB - G_LNA_dB - G_PS_dB;
P1dB_Input_LNA = P1dB_LNA_dBm - G_BPF_dB;

P1dB_Input = min([P1dB_Input_LNA, P1dB_Input_MixerI, P1dB_Input_BB1, P1dB_Input_BBA2, P1dB_Input_VGA]);

DR = P1dB_Input - MDS;

fprintf("MDS: %f dB\nP1dB_Input: %f dB\nDynamic Range: %f dB\n", MDS, P1dB_Input, DR);
