%% cargado de la señal ecg
signalStruct = load('ecg.mat');
signal = signalStruct.ecg;


[beep, beepFs] = audioread('beep1.wav');
%sound(beep, beepFs);


%% Variables útiles

Fs = 500;
L = length(signal);     % Longitud de la señal (numero de muestras)
T = L/Fs;               % Duración de la señal
t = linspace(0,T,L);


% ancho de la ventana que se va moviendo y que va integrando
% el ancho de la ventana es de 150 ms entonces se debe tener en cuenta
% la frecuencia de muestreo de la señal particular para calcular el número
% de muestras que va a contener la ventana.
n = ceil(150/((Fs.^(-1))*1000));

% mitad del número de muestras que deben haber para iniciar bien el
% algoritmo esto porque se requiere que hayan 2 latidos detectados para
% calibrar bien todos los umbrales y los estimados del intervalo RR
% esto es porque se da 360 ms para detectar un latido 

nInit =  2*ceil( (400/((Fs.^(-1))*1000)) );
RRIntervaLen = nInit/2; 

% lista de 8 intervalos RR más reciencientes 
RR1 = [];

% lista de 8 intervalos RR más recientes dentro de RR-high y RR-low
RR2 = [];

% conjunto de umbrales para la señal filtrada y la integrada
% respectivamente
thresholdSet1 = [];
thresholdSet2 = [];

% media del pico de la señal filtrada
meanPeakFromFiltSignal = 0;

% media del pico de ruido de la señal  filtrada 
meanNoisePeakFromFiltSignal = 0;

% media del pico de la señal integrada
meanPeakFromIntSignal = 0;

% media del pico del ruido de la señal integrada
meanNoisePeakFromIntSignal = 0;


% variable de lectura de la señal para simulaciones
localSignal = zeros(L,1);
localFilteredSignal = zeros(L,1);
derivative = zeros(L,1);


%% limpiar la señal    

% frecuencias de corte 
fc1 = 15; 
fc2 = 5;

% parametros para filtro pasa bajas con frecuencia de corte 15 Hz
[b1,a1] = butter(4, fc1/(Fs/2));

% parametros para filtro pasa altas con frecuencia de corte 5 Hz
[b2,a2] = butter(4, fc2/(Fs/2), 'high');


%% recibe todas las señales hasta tener cantidad suficiente

localSignal(1:n) = signal(1:n);
localFilteredSignal(1:n) = filter(b1, a1, localSignal(1:n));
localFilteredSignal(1:n) = filter(b2, a2, localFilteredSignal(1:n));

derivative(3:n-2) = customDerivative(signal(1:n), T); 
squared = derivative.^2;
integral = [zeros(n,1); movingIntegrator(squared, n)'];

%% empieza a analizar la señal 

% indices usados para calcular los intervalos RR
prevQRSindex = 0;

% indice de intervalo RR bajo, alto y miss 
RRlow = RRIntervaLen * 0.92;
RRhigh = RRIntervaLen * 1.16;
RRmiss = RRIntervaLen * 1.66;

beginning = (n+1) + 2;

for i = beginning:L-2
   
    if (i == beginning)
        prevQRSindex = i;
    end
    
    sample = signal(i);
    % filtrado de la señal para solo pasar de 5 Hz a 15 Hz 
    filteredLocalSignal = filter(b1, a1, localSignal(i-n:i));
    filteredLocalSignal = filter(b2, a2, filteredLocalSignal);
    
    localSignal(i) = sample;
    localFilteredSignal(i) = filteredLocalSignal(end);
    
    derivative(i-2) = (1/(8.*T)) .* ( (2*localFilteredSignal(i+1) + ...
                   localFilteredSignal(i+2) - localFilteredSignal(i-2)...
                   - 2*localFilteredSignal(i-1)) );
               
    squared(i-2) = derivative(i-2).^2;
    integral(i-2) = (1/n).* (sum(squared(i-2-n: i-2)));
    
    
    plt1 = subplot(3,1,1);
    plot(t(1:i-2), localSignal(1:i-2));
    title('Señal ECG Original');
    xlabel('Tiempo (s)');
    ylabel('Voltaje(mV)');

    plt2 = subplot(3,1,2);
    plot(t(1:i-2), localFilteredSignal(1:i-2));
    title('Señal ECG filtrada');
    xlabel('Tiempo (s)');
    ylabel('Voltaje(mV)');

    plt3 = subplot(3,1,3);
    plot(t(1:i-2), integral(1:i-2));
    title('Integral de la derivada al cuadrado de la señal ECG filtrada');
    xlabel('Tiempo (s)');
    ylabel('Unidades arbitrarias (u.a.)');

    if i ~= L-2
        pause(0.0000000000000001);
        cla(plt1);
        cla(plt2); 
        cla(plt3); 
    end
    

    % si se ha salido del la región RRmiss 
    if (i - prevQRSindex) > RRmiss
    
        [qrsFilteredSignal, idx1, thresholdSet1, ... 
            meanPeakFromFiltSignal, meanNoisePeakFromFiltSignal] = ... 
            detectQRScomplexSignal(localFilteredSignal(prevQRSindex:i) ...
                                , thresholdSet1 ...
                                , meanPeakFromFiltSignal...
                                , meanNoisePeakFromFiltSignal);

        [qrsIntegralSignal, idx2, thresholdSet2, ... 
            meanPeakFromIntSignal, meanNoisePeakFromIntSignal]= ... 
            detectQRScomplexSignal(integral(prevQRSindex:i) ...
            , thresholdSet2 ...
            , meanPeakFromIntSignal...
            , meanNoisePeakFromIntSignal);
        
        
        if ()
        
    end
    
    
    
    
end
















