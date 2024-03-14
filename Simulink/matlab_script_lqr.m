clc;
clear;
close all;
load('parametros4.mat');

% Equações de espaço de estados do "ramo direto"(Atuador + Dinâmica + Cinemática**)
A_estab = [-Ra/La -Kb/La 0;Ki/Jm -Bm/Jm 0;(Jr*Ki)/(Jb(1,1)*Jm) -(Jr*Bm)/(Jb(1,1)*Jm) 0];
B_estab = [1/La; 0; 0];

A_estab(3,1) = A_estab(3,1) + 1e-3;

A_apont = [-Ra/La -Kb/La 0;Ki/Jm -Bm/Jm 0;0 Jr/Jb(1,1) 0];
B_apont = [1/La;0;0];

Q_apont = [1e-4 0 0;0 5e-3 0;0 0 1e-4];
Q_estab = [1 0 0;0 10 0;0 0 20];
%Q_estab = [1 0 0;0 1 0; 0 0 1];
R_estab = 1;
R_apont = 0.1;

K_estab = lqr(A_estab,B_estab,Q_estab,R_estab);
K_apont = lqr(A_apont,B_apont,Q_apont,R_apont);

polos_apont = eig(A_apont - B_apont*K_apont);
polos_estab = eig(A_estab - B_estab*K_estab);

disp(polos_apont);
disp(polos_estab);