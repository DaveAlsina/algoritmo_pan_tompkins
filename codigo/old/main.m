
%% cargado de la señal ecg
load 'ecg.mat'

%% Variables útiles

Fs = 500;
L = 4170; % Longitud de la señal (numero de muestras)
T = L/Fs; % Duración de la señal
t = linspace(0,T,L);
f = Fs*(0:L/2)/L;
w = f*2*pi;

%% limpiar la señal    

fc1 = 15;

[b1,a1] = butter(4, fc1/(Fs/2));
filteredSignal = filter(b1,a1,ecg);

fc2 = 5;
[b2,a2] = butter(4, fc2/(Fs/2), 'high');
filteredSignal = filter(b2,a2, filteredSignal);

%% probar el moving integrator 

n= 85;

derivative = customDerivative(filteredSignal, T)';
squared = derivative.^2;
integral = movingIntegrator(squared, n)';
final = length(t)- (n+ 4); 
plot(t(1:final), integral);





