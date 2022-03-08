testRadice();
% disp(radice(-2));

function x = radice(x)
% x = radice(x)
% 
% la funziona calcola la radice di un numero 
% mediante le operazioni elementari
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
%         disp("vera: "+ vera);
%         disp("mia: "+ mia);
%         disp(abs((vera-mia)/vera));
        if (abs((vera-mia)/vera) > eps)
            disp("ITERAZIONE " + i+ " numero: " + x);
            err = err + 1; 
            disp("vera: "+ vera + "diff num rad: " + abs(x-vera^2));
            disp("mia: "+ mia + "diff num rad: " + abs(x-mia^2));
            disp ("diffe: " + abs(mia^2 - vera^2));
            disp("QUI QUALCOSA NON VA!");
        end
        x = x*10;
    end
    disp("i = " + i);
    disp("errori totali: " + err)
end