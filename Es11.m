addpath("./funcs","linsis.m");
addpath("./funcs/Es9");
[A,b] = linsis(10,1,1);
disp(mialdlt(A,b));
disp(cond(A));
[A,b] = linsis(10,10,1);
disp(cond(A));
disp(mialdlt(A,b));

% A = rand(5,5)*10;
% A = A'*A; %serve per renderla sdp
% disp(A);
% b = rand(5,1)*10;
% disp(b);
% disp("mat") 
% disp(A\b);
% disp("mia") 
% disp(mialdlt(A,b));

% A = rand(5,5)*100;
% A = A'*A;
% b = rand(5,1)*10;
% sol = mialdlt(A,b);
% writematrix(A,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(b,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')
% writematrix(sol,'matrici.csv','WriteMode','append')
% writematrix("-",'matrici.csv','WriteMode','append')