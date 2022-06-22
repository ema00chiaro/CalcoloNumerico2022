% n = 2; % n = [(1:7),9];
addpath("./funcs/Es24");

for n = 1:7
    disp("------------- "+ n +" ------------------");
    disp(weights(n));
end
disp("------------- "+ 9 + " ------------------");
disp(weights(9));
