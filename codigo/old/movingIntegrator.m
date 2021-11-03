function [integratedSignal] = movingIntegrator(signal, width)

    outputLen = length(signal); 
    integratedSignal = zeros(1, outputLen);
    
    for n = 1+width:outputLen
       
        integratedSignal(n-width) = (1/width).* (sum(signal(n-width:n)));
    end
    


end

