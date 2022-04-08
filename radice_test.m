[x,er] = testRadice(-10,10,20);
disp(x)
disp(er)
% disp(radice(-2));
% x=logspace(-10,10,20);
% for i = 1:20
%     disp("ITERAZIONE " + i);
%     disp("vera: "+ radice(x(i)));
%     disp("mia: "+ sqrt(x(i)));
% end

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

% function testRadice()
%     x=logspace(-10,10,20);
%     err = 0;
%     for i = 1:20
%         mia = radice(x(i));
%         vera = sqrt(x(i));
%         errAss = abs((vera-mia));
%         errRel = errAss/vera;
% %         disp("ITERAZIONE " + i);
% %         disp("num: " + x(i));
% %         disp ("err ass: " + errAss);
% %         disp ("err rel: " + errRel);
%         disp(x(i) +","+errRel)
%         if (errRel > eps)
%             err = err + 1; 
%             disp("vera: "+ vera + " diff num rad: " + abs(x(i)-vera^2));
%             disp("mia: "+ mia + " diff num rad: " + abs(x(i)-mia^2));
%             disp("QUI QUALCOSA NON VA!");
%         end
%     end
%     disp("i = " + i);
%     disp("errori totali: " + err)
% end


function [x,er] = testRadice(a,b,n)
% [x,er] = testRadice(a,b,n)
% 
% Funzione per testare la funzione "x = radice(x)" su n valori 
% equispaziati logaritmicamente in un certo intervallo [10^a,10^b].
% 
% INPUT:
%     a - primo esponente intervallo
%     b - secondo esponente intervallo
%     n - numero di valori equispaziati
%         logaritmicamente voluti
% 
% OUTPUT:
%     x - vettore contenente gli n valori 
%         equispaziati logaritmicamente
%     er - vettore contenente i vari errori 
%          relativi rispetto alla funzione 
%          "x = sqrt(x)" proposta da matlab

    x = logspace(a,b,n)';
    er = x;
    for i = 1:n
        er(i) = abs(radice(x(i))-sqrt((x(i))))/sqrt(x(i));
    end
    return
end