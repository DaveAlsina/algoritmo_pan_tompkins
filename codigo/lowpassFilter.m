function [filtered] = lowpassFilter(signal, n)
    
    filtered = 0;

    if (n - 12) > 0
        filtered = signal(n) - 2*signal(n-6) + signal(n-12);
    end
    
    if (n - 14) > 0
        filtered = filtered + 2*lowpassFilter(signal, n-1) - lowpassFilter(signal, n-2);
    end
    
end

