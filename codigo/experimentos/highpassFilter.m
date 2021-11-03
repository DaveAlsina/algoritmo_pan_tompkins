function [filtered] = highpassFilter(signal, n)

    filtered = 0;
    
    if (n - 32) > 0
        filtered = 32*signal(n-16) - signal(n) + signal(n-32);
    end
        
        
    if (n - 33) > 0
        filtered = filtered - highpassFilter(signal, n-1);
    end
    
end