clc
clear all
close all

load fft_data_out
% highfreq
load fft_data_out1

%%
Time = time_without;
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
xlim([0 180])
% ylim([1e-5 1e-1])
%%
Time = time_with1;
T = Time(end)-Time(end-1);                      % Sample time
Fs = 1/T;                    % Sampling frequency
dum1 = size(Time);
L = dum1(1,2);                     % Length of signal
t = (0:L-1)*T;                % Time vector
y = out_with1';
u = inp_with1';

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
U = fft(u,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1)*2*pi;
f = Fs/2*linspace(0,1,NFFT/2+1);
sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);

% Plot single-sided amplitude spectrum.
hold on
semilogy(f,2*abs(Y(1:NFFT/2+1,:)./U(1:NFFT/2+1,:)),'r') 
%%
% Time = time_with2;
% T = Time(end)-Time(end-1);                      % Sample time
% Fs = 1/T;                    % Sampling frequency
% dum1 = size(Time);
% L = dum1(1,2);                     % Length of signal
% t = (0:L-1)*T;                % Time vector
% y = out_with2';
% u = inp_with2';
% 
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% U = fft(u,NFFT)/L;
% % f = Fs/2*linspace(0,1,NFFT/2+1)*2*pi;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% sampling_time = (1-0)/(NFFT/2+1)*(Fs/2);
% 
% % Plot single-sided amplitude spectrum.
% hold on
% semilogy(f,2*abs(Y(1:NFFT/2+1,:)),'g') 
