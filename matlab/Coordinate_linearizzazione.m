clc

close all

%% Calcolo punti di equilibrio

x0 = [5,10,0]';
u0 = [0.89,0,0]';
[X, U, Y, DX] = trim('simCoordinate', x0, u0, [], 2, 1);

[A,B,C,D] = linmod('simCoordinate', X, U);

sys = ss(A,B,C,D);

%%  Funzione di trasferimento completa del sistema linearizzato (3x3)

Gc = tf(sys);

%% Funzione di trasferimento con solo gli input da disaccoppiare (2x2)

G(1,1) = Gc(1,1);
G(1,2) = Gc(1,2);
G(2,1) = Gc(2,1);
G(2,2) = Gc(2,2);

%% DELTA

DELTA(1,1) = {1};
DELTA(2,2) = {1};
DELTA(1,2) = {-G(1,2)* (1/G(1,1))};
DELTA(2,1) = {-G(2,1)* (1/G(2,2))};



