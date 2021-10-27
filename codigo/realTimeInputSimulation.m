%% cargado de la señal ecg
signalStruct = load('ecg.mat');
signal = signalStruct.ecg;

%% Variables útiles

Fs = 500;
L = 4170; % Longitud de la señal (numero de muestras)
T = L/Fs; % Duración de la señal
t = linspace(0,T,L);
f = Fs*(0:L/2)/L;
w = f*2*pi;

% ancho de la ventana que se va moviendo y que va integrando
n = 85;

%% limpiar la señal    

% frecuencias de corte 
fc1 = 15; 
fc2 = 5;

% parametros para filtro pasa bajas con frecuencia de corte 15 Hz
[b1,a1] = butter(4, fc1/(Fs/2));

% parametros para filtro pasa altas con frecuencia de corte 5 Hz
[b2,a2] = butter(4, fc2/(Fs/2), 'high');


%% simulación de leer la señal progresivamente



signaLength = length(signal); 

for i = 1:signaLength

    % el momento en que se empieza a procesar la señal es a partir de que
    % ya se tienen suficientes datos como para aplicar todas las ventanas
    
    
    if i > n
        filteredSignal = filter(b1, a1, signal);
        filteredSignal = filter(b2, a2, filteredSignal);
        
        
        derivative = customDerivative(filteredSignal, T)';
        squared = derivative.^2;
        integral = movingIntegrator(squared, n)';
        final = length(t)- (n+ 4); 
        

        plt1 = subplot(3,1,1);
        plot(t(1:i), signal(1:i));
        title('Señal ECG Original');
        xlabel('Tiempo (s)');
        ylabel('Voltaje(mV)');
        
        
        plt2 = subplot(3,1,2);
        plot(t(1:i), filteredSignal(1:i));
        title('Señal ECG filtrada');
        xlabel('Tiempo (s)');
        ylabel('Voltaje(mV)');
        
        plt3 = subplot(3,1,3);
        plot(integral(1:i));
        title('Integral de la derivada al cuadrado de la señal ECG filtrada');
        
        if i ~= signaLength
            pause(0.000000000000000000000001);
            cla(plt1);
            cla(plt2); 
            cla(plt3); 
        end
        
    end
    
    
    
    
    
end






