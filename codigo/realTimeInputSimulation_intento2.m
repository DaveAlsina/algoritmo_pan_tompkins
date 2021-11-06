%% cargado de la señal ecg
signalStruct = load('ecg2.mat');
signal = signalStruct.ecg;
Fs = signalStruct.Fs;

[beep, beepFs] = audioread('beep1.wav');
%sound(beep, beepFs);


%% Variables útiles


L = length(signal);     % Longitud de la señal (numero de muestras)
T = L/Fs;               % Duración de la señal
t = linspace(0,T,L);


% ancho de la ventana que se va moviendo y que va integrando
% el ancho de la ventana es de 150 ms entonces se debe tener en cuenta
% la frecuencia de muestreo de la señal particular para calcular el número
% de muestras que va a contener la ventana.

n = ceil(150/((Fs.^(-1))*1000));

% promedios de los 8 más grandes intervalos RR, y de los más recientes
% 8 intervalos RR, respectivamente
RRavg2 = 0;
RRavg1 = 0;

% lista de intervalos RR 
RR1 = [ ];

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



%% recibe todas las señales hasta tener cantidad suficiente

localSignal(1:n) = signal(1:n);

% para hacer el filtrado con el filtro original del paper es mejor hacerlo
% a pedazos porque si se trata hacer todo de una sola vez el algoritmo
% recursivo de filtrado tarda un siglo 

localFilteredSignal(1: floor(n/2) ) = passbandFilter(localSignal( 1:floor(n/2) ));
localFilteredSignal( floor(n/2)+1: n) = passbandFilter(localSignal( floor(n/2)+1 :n ));

derivative(3:n-2) = customDerivative(signal(1:n), T); 
squared = derivative.^2;
integral = [zeros(n,1); movingIntegrator(squared, n)'];

%% inicializa algunas variables previas que son también útiles 

% indices usados para calcular los intervalos RR
prevQRSindex = 0;
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


% variable para guardar el pico máximo que se va encontrando a lo largo del
% proceso de busqueda de picos candidatos a ser complejo QRS

maxPeakF = -Inf;
maxPeakI = -Inf;



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
    
    plt1 = subplot(4,1,1);
    plot(t(1:i-2), localSignal(1:i-2));
    title('Señal ECG Original');
    xlabel('Tiempo (s)');
    ylabel('Voltaje(mV)');

    plt2 = subplot(4,1,2);
    plot(t(1:i-2), localFilteredSignal(1:i-2));
    title('Señal ECG filtrada');
    xlabel('Tiempo (s)');
    ylabel('Voltaje(mV)');

    plt3 = subplot(4,1,3);
    plot(t(1:i-2), integral(1:i-2));
    title('Integral de la derivada al cuadrado de la señal ECG filtrada');
    xlabel('Tiempo (s)');
    ylabel('Unidades arbitrarias (u.a.)');

    plt4 = subplot(4,1,4);
    hold on
    plot(1:length(RR1), fliplr(RR1).*(Fs.^-1).*(1000), 'b--o');
    hold off
    legend('Duración estimada RR1', 'Location', 'bestoutside')
    title('Longitud RR1');
    xlabel('número de QRS detectados');
    ylabel('ancho de la ventana RR');
    
    
    if i ~= L-2
        pause(0.0000000000000001);
        cla(plt1);
        cla(plt2); 
        cla(plt3); 
        cla(plt4);
    end
    

    % Al pasar dos segundos hacemos la inicialización
    if i == twoSeconds
 
        % inicializamos el SPKI
        peaksFiltered = findpeaks(  integral(1:i) );
        [SPKI, idx] = max(peaksFiltered);
        SPKI = SPKI/3;
        
        % inicializamos el NPKI
        NPKI = mean(peaksFiltered( 1:length(peaksFiltered) ~= idx ))/2;
        
        % inicializamos el SPKF
        peaksFiltered = findpeaks( localFilteredSignal(1:i) );
        [SPKF, idx] = max(peaksFiltered);
        SPKF = SPKF/3;
        
        % inicializamos el NPKF
        NPKF = mean(peaksFiltered( 1:length(peaksFiltered) ~= idx ))/2;
               
        % Inicializamos el conjunto de umbrales para la señal filtrada 
        thresholdSetF(1) = NPKF + 0.25*(SPKF - NPKF);
        thresholdSetF(2) = thresholdSetF(1) * 0.5;
        
        % Inicializamos el conjunto de umbrales para la señal integrada
        thresholdSetI(1) = NPKI + 0.25*(SPKI - NPKI);
        thresholdSetI(2) = thresholdSetI(1) * 0.5;
        
        % Inicialización del vector de intervalos RR1
        % obseve que aquí se hace con base a la señal filtrada
        
        [peaks, peaksIdxs] = findpeaks(localFilteredSignal(1:i));
        [maxval, ~] = max(peaks); 

        % tras tener el máximo de la señal busca cuales 
        % son los picos que tienen 80% o más altura 
        % que ese pico máximo encontrado, de esta forma se busca detectar
        % los picos que pueden haber en ese intervalo inicial 

        idxs = find( peaks >= (5/10).*maxval); 
        
        % captura los indices en 'localFilteredSignal' de los picos más
        % altos
        idxs = peaksIdxs(idxs);
        
        if length(idxs) >= 2
            
            localPrevIdx = idxs(1);
            
            for i = 2:length(idxs)
                
                RR = idxs(i) - localPrevIdx;
                RR1 = [ RR RR1];
                localPrevIdx = i;
            end
            
        end
        
        prevQRSindex = idxs(end);
        
        %inicializa RRavg2
        RRavg2 = mean(RR1);
        
        RRlow = RRavg2 * 0.92;
        RRhigh = RRavg2 * 1.16;
        RRmiss = RRavg2 * 1.66;
    
    % si llevamos más de los 2 segundos y ya inicializamos entonces 
    % hacemos el proceso de detección de QRS y actualización de variables
    
    elseif i > twoSeconds
               
        isQRSFilt = false;
        isQRSInt = false;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Detecta con los umbrales si hay o no un potencial QRS %
        %   en la señal que ha sido ** FILTRADA **                %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if localFilteredSignal(i-2) > thresholdSetF(1)
        
            isQRSFilt = true;  
            SPKF = (0.125 * localFilteredSignal(i-2)) + 0.875 * SPKF;
            
            % resetea el pico máximo de la señal filtrada
            maxPeakF = -Inf;
                        
        % caso en que tiene que usar el 2do umbral para verificar 
        % si la señal pico es un QRS
        elseif (RRmiss == (i - 2 - prevQRSindex) )
            
            % si el pico máximo encontrado a lo largo de este tiempo 
            % es más grande que el segundo umbral entonces ese pico 
            % se categoriza como un complejo QRS 
            
            if maxPeakF > thresholdSetF(2)
                
                isQRSFilt = true;                
                SPKF = (0.25 * maxPeakF) + 0.75 * SPKF;
                
            % el caso en el que el pico más grande encontrado a lo largo de este
            % periodo de búsqueda, no sea más grande que el segundo umbral
            % se acutualiza el umbral de ruido
            
            else
                
                % actualizo el valor para el estimado de pico de ruido
                NPKF = (maxPeakF * 0.125) + (NPKF * 0.875); 
            end
            
            % resetea el pico máximo de la señal filtrada
            maxPeakF = -Inf;
            
        % en caso de que aún no se esté en RRmiss y el dato detectado en el
        % momento no sea mayor al umbral uno se revisa si se puede guardar
        % como máximo pico 'temporal'
        else
            
            if localFilteredSignal(i-2) > maxPeakF
                maxPeakF = localFilteredSignal(i-2);
            end
                
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Detecta con los umbrales si hay o no un potencial QRS %
        %   en la señal que ha sido ** INTEGRADA **               %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if integral(i-2) > thresholdSetI(1)
        
            isQRSInt = true;            
            SPKI = (0.125 * integral(i-2)) + 0.875 * SPKI;
            
            % resetea el pico máximo de la señal integrada
            maxPeakI = -Inf;     
            
        % caso en que tiene que usar el 2do umbral para verificar 
        % si la señal pico es un QRS
        elseif (RRmiss == (i - 2 - prevQRSindex) )
            
            if maxPeakI > thresholdSetI(2)
                
                isQRSInt = true;               
                SPKI = (0.25 * maxPeakI) + 0.75 * SPKI;
            
            else
                
                % actualizo el valor para el estimado de pico de ruido
                NPKI = (maxPeakI * 0.125) + (NPKI * 0.875); 
            end
            
            % resetea el pico máximo de la señal integrada
            maxPeakI = -Inf;
           
        % en caso de que aún no se esté en RRmiss y el dato detectado en el
        % momento no sea mayor al umbral uno se revisa si se puede guardar
        % como máximo pico 'temporal'
        
        else
            
            if integral(i-2) > maxPeakI
                maxPeakI = integral(i-2);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   DETECCIÓN DE SEÑAL, RESETEO Y ACTUALIZACIÓN DE VARS     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % setea la variable que impide que suene multiples veces 
        % la alarma en caso de detectar un QRS
        
        if ( (i - prevQRSindex) > refatoryPeriod)
            isInRefatoryPeriod = false;
        end
        
        % suena la alarma si ambos detectores dicen que hay un pico
             
        if (isQRSFilt && isQRSInt) && ( ~isInRefatoryPeriod )
            
            sound(beep, beepFs); 
            isInRefatoryPeriod = true;            
           
            % calcula la longitud del intervalo RR actual
            RR = (i - prevQRSindex);
            
            % resetea el contador
            prevQRSindex = i;
            
            % añade el RR actual a la lista
            RR1 = [RR RR1];
                       
            % saca la media de los 8 RR1 más recientes
            if length(RR1) >= 8
                RRavg1 = mean(RR1(1:8));
            else
                RRavg1 = mean(RR1( 1:length(RR1) ));
            end

            % saca la media de los 8 RR1 máximos más recientes
            if length(RR1) >= 8
                
                sortedRR1 = sort(RR1, 'descend');
                RRavg2 = mean(sortedRR1(1:8));
            else
                RRavg2 = mean(RR1( 1:length(RR1) ));
            end
           

            % indice de intervalo RR bajo, alto y miss 
            RRlow = RRavg2 * 0.92;
            RRhigh = RRavg2 * 1.16;
            RRmiss = RRavg2 * 1.66;
          
    
        end

        
    end
    


    
end
















