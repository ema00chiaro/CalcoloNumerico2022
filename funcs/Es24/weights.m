function w = weights(n)
% w = weights(n)
% 
% Funzione che restituisce i pesi di quadratura
% della formula di Newton-Cotes di grado n
% Input:
%     n - il grado della formula di Newton-Cotes
% Output:
%     w - i pesi ad essa associati
    
    if n<= 0, error("n deve avere valore positivo");end
    w = zeros(1,n+1);
    for i = 0:n
        a = poly([0:i-1,i+1:n]);
        a = [a./(n+1:-1:1), 0];
        num = polyval(a,(n));
        d = i - [0:i-1,i+1:n];
        den = prod(d);
        w(i+1) = num/den;
    end
end