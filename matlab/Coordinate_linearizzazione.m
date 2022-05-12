clc
clear
close all

%% Calcolo punti di equilibrio

x0 = [0.3491,1]';
u0 = [0.8,0,0]';
[X, U, Y, DX] = trim('simCoordinate', x0, u0, [], [], 1);

[A,B,C,D] = linmod('simCoordinate', X, U);

sys = ss(A,B,C,D);

%% 

G = tf(sys);


%OUT = G * U; 


