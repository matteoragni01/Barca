clc
clear
close all

%% Caricamento dei parametri nel workspace
parameters

%% Calcolo punti di equilibrio

y0 = 0;             
x0 = [1.179 0 0 0 y0]'; 

[X, U, Y, DX] = trim('simShip', x0, [], [], [3;5], []);

[A,B,C,D] = linmod('simShip', X, U);

sys = ss(A,B,C,D);

%% Stabilità del punto di equilibrio

autA = eig(A);
%P = lyap(A, eye(5));
J = jordan(A);

Gr = tf(sys);
% s = tf('s');
% G1 = G*(1+(1/0.05)*s);

%% Simulazione

% %open('linearizzato.slx');
% %out = sim('linearizzato.slx');
% dataNL = out.simout.Data;
% dataLIN = out.simoutLIN.Data;
% 
% uNL = dataNL(:,1);
% vNL = dataNL(:,2);
% psiNL = dataNL(:,3);
% 
% uLIN = dataLIN(:,1);
% vLIN = dataLIN(:,2);
% psiLIN = dataLIN(:,3);
% 
% 
% %% Calcolo coordinata X NL
% 
% xNL = zeros(1,length(uNL));
% xNL(1) = 0;
% for t=2:length(uNL)
%     xNL(t) = xNL(t-1) + (uNL(t)*cos(psiNL(t))-vNL(t)*sin(psiNL(t))) * 0.02;
% end
% 
% xNL = xNL';
% 
% %% Calcolo coordinata X LIN
% xLIN = zeros(1,length(uLIN));
% xLIN(1) = 0;
% for t=2:length(uLIN)
%     xLIN(t) = xLIN(t-1) + (uLIN(t)*cos(psiLIN(t))-vLIN(t)*sin(psiLIN(t))) * 0.02;
% end
% 
% xLIN = xLIN';
% 
% figure(1)
% plot(xLIN,xNL,'xb');