testRadice();

function x = radice(x)
% x = radice(x)
% 
% la funziona calcola la radice di un numero mediante le operazioni elementari
% INPUT:
%     x - il numero di cui si vuole trovare la radice
% OUTPUT:
%     x - la radice del numero
    if x < 0
        error("errore: numero negativo, impossibile fare la radice!"); 
    end
    if (x == 0 || x == 1), return, end
    sq = x;
    sq_old = 0;
    i = 0;
    while (abs(sq-sq_old) > eps)
        i = i+1;
        sq_old = sq;
        sq = (1/2)*(sq_old+(x/sq_old));
    end
    x = sq;
    return
end


function testRadice()
    inf = 1e-10;
    sup = 1e10;
    x = inf;
    err = 0;
    i = 0;
    while x <= sup
        i = i+1;
        mia = radice(x);
        vera = sqrt(x);
        disp("ITERAZIONE " + i);
        if (abs(vera - mia) > eps)
            err = err + 1; 
            disp("vera: "+ vera);
            disp("mia: "+ mia);
            disp ("diffe: " + abs(vera - mia));
            disp("QUI QUALCOSA NON VA!");
        end
        x = x*10;
    end
    disp("errori totali: " + err)
end