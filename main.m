%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
%
% Created June, 2021
% Elyes Balti
% The University of Texas at Austin
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the paper: Elyes Balti, Brian L. Evans, "Adaptive Self-Interference Cancellation for Full-Duplex
% Wireless Communication Systems," 2021. 
%
% Contact email: ebalti@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script produces the bit error rate (BER) and the spectral
% efficiency of a dual-hop full-duplex system. The BER performance are
% considered for Uncoded (without self-interference), and LDPC coded bits
% for (Interference-Free, With Interference and With Interference
% Cancellation). In the first stage, the self-interference channel is
% estimated using pilots or training frames as an input to the Steepest
% Descent method. In the second stage, the estimated self-interference
% signal is subtracted from the received uplink signal. Note that the
% uplink received signal is corrupted by the AWGN noise and the loop-back self-interference. 
% More details can be found in the paper listed above.
%% Parameters 
% SNR_dB: signal-to noise ratio (SNR) in dB
% si: self-interference power or the inverse of the signal-to-interference
% ratio (SIR)
% CHAN_LEN: number of channels taps
% CP_LEN: length of cyclic prefix
% NUM_FRAMES: number of frames (Monte Carlo iterations)
% M: modulation order
% k: number of bits per symbol
% FFT_LEN: length of FFT
% NUM_BIT: total number of bits
% FADE_VAR_1D: 1D fading variance
% kappa_dB: Rician factor in dB
% Iterations: number of iterations of the steepest descent
% alpha: parameter related to the step size
% r: coding rate
% H: parity-check matrix
% NOISE_VAR_1D: 1D VARIANCE OF AWGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc; close all
dbstop if error
SNR_dB = -10:2:30;
SI =60;
SI = db2pow(SI);
CHAN_LEN = 10;
CP_LEN = CHAN_LEN - 1;
NUM_FRAMES = 1e1; 
M = 4;
k = log2(M);
FFT_LEN = 32400; 
NUM_BIT = k*FFT_LEN;
FADE_VAR_1D = 1/2;
kappa_dB = 5; kappa = db2pow(kappa_dB);
Iterations = 1e2;
alpha = 0.125;

%% Configuring the LDPC encoder/decoder
r = 1/2;
H = dvbs2ldpc(r);
ldpcEncoder = comm.LDPCEncoder(H);% encoder
ldpcDecoder = comm.LDPCDecoder;% decoder

pskMod = comm.PSKModulator(M,'BitInput',true); % PSK modulator
pskDemod = comm.PSKDemodulator(M,'BitOutput',true,...% PSK demodulator
    'DecisionMethod','Approximate log-likelihood ratio'); 



%% Initialize the bit error vectors
PLDPCSICAN = zeros(length(SNR_dB),1);% LDPC (After SI Cancellation)
PLDPCIF = zeros(length(SNR_dB),1);% LDPC (Interference-Free)
PLDPCWI = zeros(length(SNR_dB),1);% LDPC (With Interference)

%% Initialize the spectral efficiency vectors
RIF = zeros(length(SNR_dB),1);% Interference-Free
RICAN = zeros(length(SNR_dB),1);% After SI Cancellation
RWI = zeros(length(SNR_dB),1);% With Interference


for ii=1:length(SNR_dB)
    ii
    SNR = db2pow(SNR_dB(ii));
    
    NOISE_VAR_1D = .5*2*2*FADE_VAR_1D*CHAN_LEN/(2*FFT_LEN*SNR); 
    
    %% Auxiliary vectors to evaluate the spectral efficiency for a given SNR
    C1 = zeros(NUM_FRAMES,1); C2 = zeros(NUM_FRAMES,1); C3 = zeros(NUM_FRAMES,1);  

    
    for kk=1:NUM_FRAMES
       
    %% Generate the SI channel
    hnlos = sqrt(FADE_VAR_1D)* ( randn(CHAN_LEN,1) + 1i*randn(CHAN_LEN,1)  );% NLOS
    hlos = fir1(CHAN_LEN-1,0.5)';% LOS
    hs = sqrt(SI)*( sqrt(kappa/(kappa+1))*hlos + sqrt(1/(kappa+1))*hnlos );% Aggegate SI channel

    %% Estimate the SI channel hs   
    dataIn = logical(randi([0 1],FFT_LEN,1));
    encData = ldpcEncoder(dataIn); 
    txSig = pskMod(encData);
    TX_SIG = ifft(txSig);
    CP = TX_SIG(end-CP_LEN+1:end);
    PILOT = [CP; TX_SIG];
    hest = SteepestDescent(PILOT,hs,Iterations,alpha,NOISE_VAR_1D);   
    
    %% Transmit Data Frame (Downlink data leaked to loop-back SI path)
    dataDownlink = logical(randi([0 1],FFT_LEN,1));
    encData = ldpcEncoder(dataDownlink); 
    txSig = pskMod(encData);
    TX_SIG = ifft(txSig);
    CP = TX_SIG(end-CP_LEN+1:end);
    TxFrameDownlink = [CP; TX_SIG];  
        
    %% Uplink Transmit Data Frame
    data = logical(randi([0 1],FFT_LEN,1));
    encData = ldpcEncoder(data); 
    txSig = pskMod(encData);
    TX_SIG = ifft(txSig);
    CP = TX_SIG(end-CP_LEN+1:end);
    TxFrame = [CP; TX_SIG];
    
    %% Generate the uplink channel and the noise
    hu = sqrt(FADE_VAR_1D)* ( randn(CHAN_LEN,1) + 1i*randn(CHAN_LEN,1)  );
    noise = randn(FFT_LEN + CP_LEN+CHAN_LEN-1,1) + 1i*randn(FFT_LEN + CP_LEN+CHAN_LEN-1,1);
    
   %% LDPC Coded Reception (With Interference)  
     yu = conv(hu,TxFrame) + conv(hs,TxFrameDownlink) + sqrt(NOISE_VAR_1D) *noise;
     yu(1:CP_LEN) = [];
     yu = yu(1:FFT_LEN);
     Yu = fft(yu,FFT_LEN);
     Hu = fft(hu,FFT_LEN);
     RX_EQ = Yu./Hu;
     demodSig= pskDemod(RX_EQ);
     rxBits = ldpcDecoder(demodSig);
     nErrors = biterr(rxBits,data);
     PLDPCWI(ii) = PLDPCWI(ii) + nErrors;
    
    %% LDPC Coded Reception (With SI Cancellation)
    yu = conv(hu,TxFrame) + conv(hs,TxFrameDownlink) - conv(hest,TxFrameDownlink) + sqrt(NOISE_VAR_1D) *noise;
    yu(1:CP_LEN) = [];
    yu = yu(1:FFT_LEN);
    Yu = fft(yu,FFT_LEN);
    RX_EQ = Yu./Hu;
    demodSig= pskDemod(RX_EQ);
    rxBits = ldpcDecoder(demodSig);
    nErrors = biterr(rxBits,data);
    PLDPCSICAN(ii) = PLDPCSICAN(ii) + nErrors;
    
    %% LDPC Coded Reception (Without Interference)
    yu = conv(hu,TxFrame) + sqrt(NOISE_VAR_1D) *noise;
    yu(1:CP_LEN) = [];
    yu = yu(1:FFT_LEN);
    Yu = fft(yu,FFT_LEN);
    RX_EQ = Yu./Hu;
    demodSig= pskDemod(RX_EQ);
    rxBits = ldpcDecoder(demodSig);
    nErrors = biterr(rxBits,data);
    PLDPCIF(ii) = PLDPCIF(ii) + nErrors;
    
        
     %% Spectral Efficiency
     
     % Apply the FFT to the SI channels
      Hs = fft(hs,FFT_LEN); 
      Hest = fft(hest,FFT_LEN);
      EE = Hs - Hest;
      
     for nn=1:FFT_LEN
     C1(kk) = C1(kk) + log2(1 + abs(Hu(nn)).^2/FFT_LEN/ (2*NOISE_VAR_1D) );
     C2(kk) = C2(kk) + log2(1 + abs(Hu(nn)).^2/FFT_LEN/ (2*NOISE_VAR_1D + abs(EE(nn)).^2/FFT_LEN ));
     C3(kk) = C3(kk) + log2(1 + abs(Hu(nn)).^2/FFT_LEN/ (2*NOISE_VAR_1D + abs(Hs(nn)).^2/FFT_LEN  ) );
     end
     
     end
  RIF(ii) = mean(C1)/FFT_LEN;
  RICAN(ii) = mean(C2)/FFT_LEN;
  RWI(ii) = mean(C3)/FFT_LEN;

end


%% Averaging over the total number of frames and bits
PLDPCWI = PLDPCWI/(NUM_BIT*NUM_FRAMES);
PLDPCSICAN = PLDPCSICAN/(NUM_BIT*NUM_FRAMES);
PLDPCIF = PLDPCIF/(NUM_BIT*NUM_FRAMES);

%% BER (Uncoded Without Interference)
PUNCIF = berfading(SNR_dB,'qam',M,1);



 figure; hold on;
 set(gca,'yscale','log')
 plot(SNR_dB,PUNCIF,'linewidth',2)
 plot(SNR_dB,PLDPCSICAN,'linewidth',2)
 plot(SNR_dB,PLDPCIF,'linewidth',2)
 plot(SNR_dB,PLDPCWI,'linewidth',2)
 xlabel('SNR (dB)')
 ylabel('BER')
 legend('Uncoded (Interference-Free)','LDPC (After SI Cancellation)','LDPC (Interference-Free)','LDPC (With Interference)','location','southwest')


 figure; hold on;
 plot(SNR_dB,RIF,'linewidth',2)
 plot(SNR_dB,RICAN,'linewidth',2)
 plot(SNR_dB,RWI,'linewidth',2)
 xlabel('SNR (dB)')
 ylabel('Spectral Efficiency (bits/s/Hz)')
 legend('Interference-Free','After SI Cancellation','With Interference','location','northwest')



