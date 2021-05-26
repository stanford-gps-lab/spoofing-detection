function R = Rcorr(tau)
    %R = Rcorr(tau)
    %   Calculates the value of the triangle autocorrelation
    %   function at a value of tau chips from the prompt.

    R = max(1 - abs(tau), 0);

end