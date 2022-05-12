clc
clear
close all

%% Caricamento dei parametri nel workspace
parameters

%% Simulazione
out = sim("simTurning_35_Step");

%% Output
data = out.simout.Data;

psi = data(:,1);
xaxis = data(:,2);
yaxis = data(:,3);

%% Grafico
figure('Name', "Ship-NONLIN")
xprimo = xaxis/p.Lpp;
yprimo = yaxis/p.Lpp;

plot(yprimo,xprimo)
title("Turning manoeuvre (35° - Step)")
xlabel("Position - X");
ylabel("Position - Y");





