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

RRavg2 = ceil( (400/((Fs.^(-1))*1000)) );
RRavg1 = 0;

% lista de 8 intervalos RR más reciencientes 
RR1 = zeros(8,1);

% lista de 8 intervalos RR más recientes dentro de RR-high y RR-low
RR2 = [ RRavg2 ];

% conjunto de umbrales para la señal filtrada y la integrada
% respectivamente
thresholdSetF = [];
thresholdSetI = [];

% media del pico de la señal filtrada
SPKF = 0;

% media del pico de ruido de la señal  filtrada 
NPKF = 0;

% media del pico de la señal integrada
SPKI = 0;

% media del pico del ruido de la señal integrada
NPKI = 0;


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
%localFilteredSignal(1:n) = filter(b1, a1, localSignal(1:n));
%localFilteredSignal(1:n) = filter(b2, a2, localFilteredSignal(1:n));
localFilteredSignal(1: floor(n/2) ) = passbandFilter(localSignal( 1:floor(n/2) ));
localFilteredSignal( floor(n/2)+1: n) = passbandFilter(localSignal( floor(n/2)+1 :n ));

derivative(3:n-2) = customDerivative(signal(1:n), T); 
squared = derivative.^2;
integral = [zeros(n,1); movingIntegrator(squared, n)'];

%% inicializa algunas variables previas que son también útiles 

% indices usados para calcular los intervalos RR
prevQRSindex = 0;

% indice de intervalo RR bajo, alto y miss 
RRlow = RRavg2 * 0.92;
RRhigh = RRavg2 * 1.16;
RRmiss = RRavg2 * 1.66;

beginning = (n+1) + 2;

% tiempo que debe esperar para inicializar la detección de umbrales
twoSeconds = ceil(2000/((Fs.^(-1))*1000));

% tiempo en el que tras ocurrir un latido no puede ocurrir otro 
% (lo dice el paper, son 200 ms, en este caso se pasan esos segundos 
% a número de muestras en que se alcanzan esos segundos). 

refatoryPeriod = ceil(200/((Fs.^(-1))*1000));

% variable para indicar que está en periodo de espera para poder detectar 
% otro QRS

isInRefatoryPeriod = false;


%% empieza a analizar la señal en tiempo real

for i = beginning:L-2
   
    % inicializa prevQRSindex
    if (i == beginning)
        prevQRSindex = i;
    end
    
    %obtiene la muestra actual
    sample = signal(i);
    localSignal(i) = sample;
    
    % filtrado de la señal para solo pasar de 5 Hz a 15 Hz 
    % con passbandFilter se implementa literalmente el filtro 
    % de pam-tompkins en su paper 
    
    filteredLocalSignal = passbandFilter(localSignal(i-34:i));
    localFilteredSignal(i) = filteredLocalSignal(end);
    
    derivative(i-2) = (1/(8.*T)) .* ( (2*localFilteredSignal(i+1) + ...
                   localFilteredSignal(i+2) - localFilteredSignal(i-2)...
                   - 2*localFilteredSignal(i-1)) );
               
    squared(i-2) = derivative(i-2).^2;
    integral(i-2) = (1/n).* (sum(squared(i-2-n: i-2)));
    
    % Hace los plots (observamos que el cuello de botella está en la 
    % graficación no en el procesamiento de los datos )
    
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
    

    % Al pasar dos segundos hacemos la inicialización
    if i == twoSeconds
 
        % inicializamos el SPKF
        peaksFiltered = findpeaks( localFilteredSignal(1:i) );
        [SPKF, idx] = max(peaksFiltered);
        SPKF = SPKF/3;
        
        % inicializamos el NPKF
        NPKF = mean(peaksFiltered( 1:length(peaksFiltered) ~= idx ))/2;
        
        % inicializamos el SPKI
        peaksFiltered = findpeaks(  integral(1:i) );
        [SPKI, idx] = max(peaksFiltered);
        SPKI = SPKI/3;
        
        % inicializamos el NPKI
        NPKI = mean(peaksFiltered( 1:length(peaksFiltered) ~= idx ))/2;
               
        % Inicializamos el conjunto de umbrales para la señal filtrada 
        thresholdSetF(1) = NPKF + 0.25*(SPKF - NPKF);
        thresholdSetF(2) = thresholdSetF(1) * 0.5;
        
        % Inicializamos el conjunto de umbrales para la señal integrada
        thresholdSetI(1) = NPKI + 0.25*(SPKI - NPKI);
        thresholdSetI(2) = thresholdSetI(1) * 0.5;
        
        prevQRSindex = idx;
        sound(beep, beepFs);
    
    % si llevamos más de los 2 segundos y ya inicializamos entonces 
    % hacemos el proceso de detección de QRS y actualización de variables
    
    elseif i > twoSeconds
        
        isQRSFilt = false;
        isQRSInt = false;
        
        %   Detecta con los umbrales si hay o no un potencial QRS 
        %   en la señal que ha sido filtrada
        
        if localFilteredSignal(i-2) > thresholdSetF(1)
        
            isQRSFilt = true;
            % añade el nuevo intervalo RR a la lista RR1
            RR1 = [(i - prevQRSindex); RR1];
            
            % añade el nuevo intervalo RR a la lista RR2
            if (RR1(1) >= RRlow ) && (RR1(1) <= RRhigh )
                RR2 = [RR1(1); RR2];
            end
            
            prevQRSindex = i-2;
            SPKF = (0.125 * localFilteredSignal(i-2)) + 0.875 * SPKF;
                        
        % caso en que tiene que usar el 2do umbral para verificar 
        % si la señal pico es un QRS
        elseif (RRmiss == (i - 2 - prevQRSindex) )
            
            % debería entonces regresarme en el pasado a revisar 
            % si hay una señal que cumpla, eso es paila, toca mejorar el
            % código
            if localFilteredSignal(i - 2) > thresholdSetF(2)
                
                isQRSFilt = true;
                %añade el nuevo intervalo RR a la lista 
                RR1 = [(i - prevQRSindex); RR1];
                
                %añade el nuevo intervalo RR a la lista RR2
                if (RR1(1) >= RRlow ) && (RR1(1) <= RRhigh )
                    RR2 = [RR1(1); RR2];
                end

                prevQRSindex = i-2; 
                SPKF = (0.25 * localFilteredSignal(i-2)) + 0.75 * SPKF;
                
            end
            
        end
        

        %   Detecta con los umbrales si hay o no un potencial QRS 
        %   en la señal que ha sido integrada
        
        if integral(i-2) > thresholdSetI(1)
        
            isQRSInt = true;
            prevQRSindex = i-2;
            
            SPKI = (0.125 * integral(i-2)) + 0.875 * SPKI;
            
            
        % caso en que tiene que usar el 2do umbral para verificar 
        % si la señal pico es un QRS
        elseif (RRmiss == (i - 2 - prevQRSindex) )
            
            if integral(i - 2) > thresholdSetI(2)
                isQRSInt = true;
                prevQRSindex = i-2;
                
                SPKI = (0.25 * integral(i-2)) + 0.75 * SPKI;
            end
            
        end
        
        
        % setea la variable que impide que suene multiples veces 
        % la alarma en caso de detectar un QRS
        
        if ( (i - prevQRSindex) > refatoryPeriod)
            isInRefatoryPeriod = false;
        end
        
        % suena la alarma si ambos detectores dicen que hay un pico
             
        if (isQRSFilt && isQRSInt) && ( ~isInRefatoryPeriod )
            
           sound(beep, beepFs); 
           isInRefatoryPeriod = true;
        end
        
        % saca la media de los 8 RR1 más recientes
        if length(RR1) >= 8
            RRavg1 = mean(RR1(1:8));
        else
            RRavg1 = mean(RR1( 1:length(RR1) ));
        end
        
        % saca la media de los 8 RR2 más recientes
        if length(RR2) >= 8
            RRavg2 = mean(RR2(1:8));
        else
            RRavg2 = mean(RR2( 1:length(RR2) ));
        end

        % indice de intervalo RR bajo, alto y miss 
        RRlow = RRavg2 * 0.92;
        RRhigh = RRavg2 * 1.16;
        RRmiss = RRavg2 * 1.66;
        
    end
    


    
end
















