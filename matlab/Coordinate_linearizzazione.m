clc
clear
close all

%% Calcolo punti di equilibrio

x0 = [0,0,0.3491]';

[X, U, Y, DX] = trim('simCoordinate', x0, [], [], [], []);

[A,B,C,D] = linmod('simCoordinate', X, U);

sys = ss(A,B,C,D);

%% 

G = tf(sys);

% Pagina 730 libro pdf
% Pagina 433 libro cartaceo
%OUT = G * U; 


