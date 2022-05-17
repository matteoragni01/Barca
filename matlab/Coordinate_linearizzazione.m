clc
clear
close all

%% Calcolo punti di equilibrio

x0 = [0,0,0.3491]';

[X, U, Y, DX] = trim('simCoordinate', x0, [], [], [], []);

[A,B,C,D] = linmod('simCoordinate', X, U);

sys = ss(A,B,C,D);

%%  Funzione di trasferimento completa del sistema linearizzato (3x3)

Gc = tf(sys);

%% Funzione di trasferimento con solo gli input da disaccoppiare (2x2)

G(1,1) = Gc(1,1);
G(1,2) = Gc(1,2);
G(2,1) = Gc(2,1);
G(2,2) = Gc(2,2);

%% Funzione di trasferimento in cui gli input sono disaccopiati

Gd(1,1) = G(1,1);
Gd(1,2) = 0;
Gd(2,1) = 0;
Gd(2,2) = G(2,2);

%% DELTA

D= inv(G)*Gd;


