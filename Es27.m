addpath("./funcs/Es24");
addpath("./funcs/Es26");
tols = [10^-2;10^-3;10^-4;10^-5;10^-6];
ris = zeros(1,5);
vals = zeros(1,5);
errs = zeros(1,5);
% for i = 1:9
%     if i == 8
%     else
%         disp("ITERAZIONE-------- " + i);
%         [ris(1),errs(1),vals(1)] = composita(@functionToPass,0,1,i,10^-2);
%         [ris(2),errs(2),vals(2)] = composita(@functionToPass,0,1,i,10^-3);
%         [ris(3),errs(3),vals(3)] = composita(@functionToPass,0,1,i,10^-4);
%         [ris(4),errs(4),vals(4)] = composita(@functionToPass,0,1,i,10^-5);
%         [ris(5),errs(5),vals(5)] = composita(@functionToPass,0,1,i,10^-6);
%         disp(table(tols,ris',errs',vals'));
%     end
% end
for i = 1:9
    if i == 8
    else
        disp("ITERAZIONE-------- " + i);
        [ris(1),errs(1),vals(1)] = compositaNew(@functionToPass,0,1,i,10^-2);
        [ris(2),errs(2),vals(2)] = compositaNew(@functionToPass,0,1,i,10^-3);
        [ris(3),errs(3),vals(3)] = compositaNew(@functionToPass,0,1,i,10^-4);
        [ris(4),errs(4),vals(4)] = compositaNew(@functionToPass,0,1,i,10^-5);
        [ris(5),errs(5),vals(5)] = compositaNew(@functionToPass,0,1,i,10^-6);
        disp(table(tols,ris',errs',vals'));
    end
end

% composita(@functionToPass,0,1,2,10^-4);

function y =  functionToPass(x)
    y = sin(1./(0.1+x));
%     y = exp(3.*x);
end