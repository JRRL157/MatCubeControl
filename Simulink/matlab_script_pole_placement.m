clc;
clear;
close all;
load('parametros4.mat');

% Equações de espaço de estados do "ramo direto"(Atuador + Dinâmica + Cinemática**)
A_estab = [-Ra/La -Kb/La 0;Ki/Jm -Bm/Jm 0;(Jr*Ki)/(Jb(1,1)*Jm) -(Jr*Bm)/(Jb(1,1)*Jm) 0];
B_estab = [1/La; 0; 0];

A_apont = [-Ra/La -Kb/La 0;Ki/Jm -Bm/Jm 0;0 Jr/Jb(1,1) 0];
B_apont = [1/La;0;0];

%Calculando os polos considerando a forma padrão de segunda ordem e os
%parâmetros ts(Acomodação), Máximo sobressinal(Mp) e 0 < p < 1(variação do
%tempo de acomodação)
%wn_apont,sigma_apont,wd_apont,csi_apont = poles(Mp_apont,p_apont,ts_apont);
%wn_estab,sigma_estab,wd_estab,csi_estab = poles(Mp_estab,p_estab,ts_estab);
ret_apont = poles(Mp_apont,p_apont,ts_apont);
ret_estab = poles(Mp_estab,p_estab,ts_estab);

wn_apont = ret_apont(1);
sigma_apont = ret_apont(2);
wd_apont = ret_apont(3);
csi_apont = ret_apont(4);

wn_estab = ret_estab(1);
sigma_estab = ret_estab(2);
wd_estab = ret_estab(3);
csi_estab = ret_estab(4);

if sigma_apont < 0
    sigma_apont = -sigma_apont;
end

if wd_apont < 0
    wd_apont = -wd_apont;
end

if sigma_estab < 0
    sigma_estab = -sigma_estab;
end

if wd_estab < 0
    wd_estab = -wd_estab;
end


%É necessário adicionar mais polos no infinito(Não dominantes) para ficar
%com a mesma dimensão n da matriz A

polos_apont = [-sigma_apont+1i*wd_apont -sigma_apont-1i*wd_apont -sigma_apont*10];
polos_estab = [-sigma_estab+1i*wd_estab -sigma_estab-1i*wd_estab -sigma_estab*0.01];

K_apont = place(A_apont,B_apont,polos_apont);
K_estab = place(A_estab,B_estab,polos_estab);

placed_poles_estab = eig(A_estab -B_estab*K_estab);
placed_poles_apont = eig(A_apont - B_apont*K_apont);

function y = poles(Mp,p,ts)
    csi = sqrt(log(Mp)*log(Mp)/(pi^2 + log(Mp)^2));
    wn = (log(100) - log(100*p) - log(sqrt(1 - csi^2)))/ts*csi;
    sigma = wn*csi;
    wd = wn*sqrt(1-csi^2);
    y = [wn,sigma,wd,csi];
end