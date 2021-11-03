%% cargado de la señal ecg
signalStruct = load('ecg.mat');
signal = signalStruct.ecg;


[beep, beepFs] = audioread('beep1.wav');
%sound(beep, beepFs);

%% Variables útiles

Fs = 500;
L = length(signal); % Longitud de la señal (numero de muestras)
T = L/Fs; % Duración de la señal
t = linspace(0,T,L);
f = Fs*(0:L/2)/L;
w = f*2*pi;

% ancho de la ventana que se va moviendo y que va integrando
% el ancho de la ventana es de 150 ms entonces se debe tener en cuenta
% la frecuencia de muestreo de la señal particular para calcular el número
% de muestras que va a contener la ventana.
n = ceil(150/((Fs.^(-1))*1000));

% mitad del número de muestras que deben haber para iniciar bien el
% algoritmo esto porque se requiere que hayan 2 latidos detectados para
% calibrar bien todos los umbrales y los estimados del intervalo RR
% esto es porque se da 360 ms para detectar un latido 

nInit =  ceil( (400/((Fs.^(-1))*1000)) );
RRIntervaLen = nInit *2; 

%% limpiar la señal    

% frecuencias de corte 
fc1 = 15; 
fc2 = 5;

% parametros para filtro pasa bajas con frecuencia de corte 15 Hz
[b1,a1] = butter(4, fc1/(Fs/2));

% parametros para filtro pasa altas con frecuencia de corte 5 Hz
[b2,a2] = butter(4, fc2/(Fs/2), 'high');


%% variables para la simulación

signaLength = length(signal); 

% filtrado de la señal para solo pasar de 5 Hz a 15 Hz 
filteredSignal = filter(b1, a1, signal);
filteredSignal = filter(b2, a2, filteredSignal);

% derivada de la señal (note como se añaden unos ceros al principio de la
% señal para compensar el número de datos y que la derivada quede de igual
% largo a la señal filtrada)

derivative = [zeros(4,1); customDerivative(filteredSignal, T)'];

% cuadrado de la señal derivada
squared = derivative.^2;

% calculo de la integral del cuadrado de la derivada de la señal filtrada
integral = [zeros(n,1); movingIntegrator(squared, n)'];

% conjunto de umbrales para la señal filtrada
% contiene el umbral 1 y 2 respectivamente
thresholdSet1 = [];

% conjunto de umbrales para la señal integrada
% contiene el umbral 1 y 2 respectivamente
thresholdSet2 = [];

meanPeakFromFiltSignal = 0;
meanNoisePeakFromFiltSignal = 0;

meanPeakFromIntSignal = 0;
meanNoisePeakFromIntSignal = 0;

result1 = false;
result2 = false;

%% simulación de leer la señal progresivamente

for i = 1:signaLength

    % el momento en que se empieza a procesar la señal es a partir de que
    % ya se tienen suficientes datos como para aplicar todas las ventanas

    if i > n
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
        plot(t(1:i), integral(1:i));
        title('Integral de la derivada al cuadrado de la señal ECG filtrada');
        xlabel('Tiempo (s)');
        ylabel('Unidades arbitrarias (u.a.)');
        
        if i ~= signaLength
            pause(0.0000000000000001);
            cla(plt1);
            cla(plt2); 
            cla(plt3); 
        end
        
    end
    
    %% 
    if ( mod(i-1, RRIntervaLen) == 0) && (i-1 > 1)
    
    [result1, idx1, thresholdSet1, meanPeakFromFiltSignal, meanNoisePeakFromFiltSignal] = ...
            detectQRScomplexSignal(filteredSignal(i-RRIntervaLen:i), thresholdSet1, ... 
            meanPeakFromFiltSignal, meanNoisePeakFromFiltSignal);

    [result2, idx2, thresholdSet2, meanPeakFromIntSignal, meanNoisePeakFromIntSignal] = ...
                detectQRScomplexSignal(integral(i-RRIntervaLen:i), thresholdSet2, ... 
                meanPeakFromIntSignal, meanNoisePeakFromIntSignal);
    
    end
    
    if (result2 && result1) && (i > nInit)
                              
          sound(beep, beepFs);
          
    end
    
    result1 = false;
    result2 = false;
end






