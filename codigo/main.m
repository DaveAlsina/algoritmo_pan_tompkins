
%% cargado de la señal ecg
load 'ecg.mat'

ecg_fft = fft(ecg);

Fs = 500;
L = 4170; % Longitud de la señal (numero de muestras)
T = L/Fs; % Duración de la señal
t = linspace(0,T,L);
f = Fs*(0:L/2)/L;
w = f*2*pi;

ecg_fft2 = abs(ecg_fft/L);
ecg_fft3 = ecg_fft2(1:L/2+1);
ecg_fft3(2:end-1) = 2*ecg_fft3(2:end-1);

plot(f, ecg_fft3)

%% limpiar la señal    

ecg_fft3( (f < 5) | (f >15) )  = 0
plot(f, ecg_fft3)

