function [outputArg1,outputArg2] = movingIntegrator(signal, width)

    outputLen = length(signal) - width; 
    
    for n = 1:outputLen
       
        integratedSignal(n) = sum(signal(n:n+width));
    end
    


end

