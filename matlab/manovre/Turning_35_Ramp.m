clc
clear
close all

%% Caricamento dei parametri nel workspace
parameters

%% Simulazione
out = sim("simTurning_35_Ramp");

%% Output
data = out.simout.Data;

psi = data(:,1);
x = data(:,2);
y = data(:,3);

%% Grafico
figure('Name', "Ship-NONLIN")
xprimo = x/p.Lpp;
yprimo = y/p.Lpp;

plot(yprimo,xprimo)
title("Turning manoeuvre (35° - Ramp)")
xlabel("Position - X");
ylabel("Position - Y");


