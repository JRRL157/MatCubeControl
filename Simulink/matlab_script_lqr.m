clc;
clear;
close all;
load('parametros3.mat');

% Equações de espaço de estados do "ramo direto"(Atuador + Dinâmica + Cinemática**)
A_estab = [-Ra/La -Kb/La 0;Ki/Jm -Bm/Jm 0;(Jr*Ki)/(Jb(1,1)*Jm) -(Jr*Bm)/(Jb(1,1)*Jm) 0];
B_estab = [1/La; 0; 0];

A_apont = [-Ra/La -Kb/La 0;Ki/Jm -Bm/Jm 0;0 Jr/Jb(1,1) 0];
B_apont = [1/La;0;0];

K_estab = lqr(A_estab,B_estab,Q_estab,R_estab);
K_apont = lqr(A_apont,B_apont,Q_apont,R_apont);

polos_apont = eig(A_apont - B_apont*K_apont);
polos_estab = eig(A_estab - B_estab*K_estab);

disp(polos_apont);
disp(polos_estab);