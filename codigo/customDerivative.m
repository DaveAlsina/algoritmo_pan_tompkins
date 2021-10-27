function [derivative] = customDerivative(signal, period)

    n = length(signal);
    derivative = zeros(1, n-4);

    for i = 3:n-2
        derivative(i-2) = (1/(8.*period)) .* ( (2*signal(i+1) + signal(i+2) - signal(i-2) - 2*signal(i-1)) );  
    end

end

