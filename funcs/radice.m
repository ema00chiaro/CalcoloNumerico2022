function x = radice(x)
% x = radice(x)
% 
% la funziona calcola la radice di un numero 
% mediante le operazioni algebriche elementari
% INPUT:
%     x - il numero di cui si vuole trovare la radice
% OUTPUT:
%     x - la radice del numero
    if x < 0
        error("errore: numero negativo, " + ...
            "impossibile calcolare la radice!"); 
    end
    if (x == 0 || x == 1), return, end
    
    sq_old = x;
    sq = x/2 + 1/2; % prima iterazione
    while ( abs(sq - sq_old) > eps)
        sq_old = sq;
        sq = (1/2)*(sq_old+(x/sq_old));
    end
    x = sq;
    return
end