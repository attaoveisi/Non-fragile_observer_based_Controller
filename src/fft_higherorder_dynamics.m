clc
clear all
close all

load fft_data_higher_order_dynamics


%%
Time = time;
T = Time(end)-Time(end-1);                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(Time);
L = dum1(1,2);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = out_without';
u = inp_without';

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1)*2*pi;
f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

% Plot single-sided amplitude spectrum.
% figure;
semilogy(f,2*abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:))) 
% hold on
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (rad/sec)')
% ylabel('|Y(f)|')
% xlim([0 35*2*pi])
xlim([180 240])
% ylim([1e-5 1e-1])
