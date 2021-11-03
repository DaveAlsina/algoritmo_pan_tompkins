function [qrsIntegralSignal, idx, thresholdSet1, meanSignalPeak, meanNoisePeak] = detectQRScomplexSignal(signal, thresholdSet1, meanSignalPeak, meanNoisePeak)

    
    n = length(signal);
    [maxval, idx] = max(signal);
    
    % variables para detectar si es o no un complejo QRS
    qrsIntegralSignal = false;
    
    % si los umbrales están vacíos quiere decir que el algoritmo hasta 
    % ahora está empezando por lo que debe tomar uno de los valores 
    % máximos en el RR como pico
    
    if isempty(thresholdSet1)
        
        meanSignalPeak = 0.875*meanSignalPeak + 0.125*maxval;
    
    % si el valor maximo en este intervalo RR es mayor que el umbral 1 
    % entonces debe ser un QRS
    
    elseif maxval > thresholdSet1(1) 
        
       qrsIntegralSignal = true;
       meanSignalPeak = 0.875*meanSignalPeak + 0.125*maxval;
       
    % de lo contrario hay que revisar si este latido tiene un valor más
    % alto que el del segundo umbral, de ser así se reajusta la
    % sensibilidad del umbral superior y se acepta ese valor como complejo
    % QRS
    
    elseif maxval > thresholdSet1(2) 
            
       qrsIntegralSignal = true;
       meanSignalPeak = 0.75*meanSignalPeak + 0.25*maxval;

    % en caso de no ser por ningún medio un candidato a complejo QRS 
    % se toma este valor máximo como si fuese ruido y con él se actualizan
    % los umbrales de ruido
    else
        
        meanNoisePeak = 0.875*meanNoisePeak + 0.125*maxval;
        
    end
            
    % después de todo lo anterior se actualizan los umbrales 
    
    thresholdSet1(1) = meanNoisePeak + 0.25*(meanSignalPeak - meanNoisePeak);
    thresholdSet1(2) = thresholdSet1(1)*0.5;
    
end

