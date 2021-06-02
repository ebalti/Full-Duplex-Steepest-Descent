%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SteepestDescent.m
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
% This script implements the steepest descent method applied to estimate
% the self-interference channel in the fist stage (training cycle). The
% implementation is based on the steps described in the book: 
% S. Haykin, "Adaptive Filter Theory (5th Ed.)," 2013-2014.
%% Parameters 
% PILOT: pilot or training frame used for channel estimation
% hs: self-interference channel
% Iterations: number of iterations of the steepest descent
% alpha: parameter related to the step size
% NOISE_VAR_1D: 1D VARIANCE OF AWGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = SteepestDescent(PILOT,hs,Iterations,alpha,NOISE_VAR_1D)
    CHAN_LEN = length(hs);
    FFT_LEN = length(PILOT) - CHAN_LEN +1;
    CP_LEN = CHAN_LEN-1;
    
    %% Construct the auto-correlation matrix of the pilots
    AUTO_CORR_VEC = xcorr(PILOT,PILOT,'unbiased').';
    MID_POINT = (length(AUTO_CORR_VEC)+1)/2;
    C = AUTO_CORR_VEC(MID_POINT:MID_POINT+CHAN_LEN-1); % FIRST COLUMN OF TOEPLITZ MATRIX
    R = fliplr(AUTO_CORR_VEC(MID_POINT-CHAN_LEN+1:MID_POINT)); % FIRST ROW OF TOEPLITZ MATRIX
    Rvv_MATRIX = toeplitz(C,R);

    MAX_STEP_SIZE = 1/real(max(eig(Rvv_MATRIX)));
    mu= alpha*MAX_STEP_SIZE;  
    
    %% Generate the noise and the SI signal
    noise = randn(FFT_LEN + CP_LEN+CHAN_LEN-1,1) + 1i*randn(FFT_LEN + CP_LEN+CHAN_LEN-1,1);
    Xsi = conv(hs,PILOT) + sqrt(NOISE_VAR_1D) *noise;
    
    %% Construct the cross-correlation matrix between the pilots and the noisy SI signal
    CROSS_CORR_VEC = xcorr(Xsi(1:length(PILOT)),PILOT,'unbiased');
    MID_POINT = (length(CROSS_CORR_VEC)+1)/2;
    IP_CROSS_CORR_VEC = CROSS_CORR_VEC(MID_POINT:MID_POINT+CHAN_LEN-1);
    
    %% Initialization of the estimated SI channel
    hest = zeros(CHAN_LEN,1);

    %% Steepest Descent Adaptation Cycle
    for i1 = 1:Iterations
     hest = hest + mu*(IP_CROSS_CORR_VEC - Rvv_MATRIX*hest);  
    end

    output=hest;
end