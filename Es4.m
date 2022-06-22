addpath("./funcs/Es4");

[x,er] = testRadice(-10,10,20);
disp(er);

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