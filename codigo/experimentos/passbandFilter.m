function [filteredSignal] = passbandFilter(signal)

    n = length(signal);
    filteredSignal = zeros(n,1);
    
    for i = 1:n
        
        filteredSignal(i) = lowpassFilter(signal, i);
        filteredSignal(i) = lowpassFilter(filteredSignal, i);
        
        disp("iteracion: ");
        disp(i);
    end


end

