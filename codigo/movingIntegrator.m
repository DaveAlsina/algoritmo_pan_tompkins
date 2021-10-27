function [integratedSignal] = movingIntegrator(signal, width)

    outputLen = length(signal) - width; 
    integratedSignal = zeros(1, outputLen);
    
    for n = 1:outputLen
       
        integratedSignal(n) = (1/width).* (sum(signal(n:n+width)));
    end
    


end

