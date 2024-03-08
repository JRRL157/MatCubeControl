clc;
clear;
close all;
load('parametros4.mat');

% Funções de transferência do ramo direto A(S)*P(s)
K_tf_forward_apont = (Ki*Jb(1,1))/(Jm*La*Jr);
K_tf_forward_estab = (Ki*Jb(1,1))/(Jm*La*Jr);

TF_forward_estab = tf([K_tf_forward_apont],[1 (Ra/La + Bm/Jm) (Ki*Kb + Bm*Ra)/(Jm*La)]);
TF_forward_apont = tf([K_tf_forward_estab],[1 (Ra/La + Bm/Jm) (Ki*Kb + Bm*Ra)/(Jm*La) 0]);
polos_tf_apont = pole(TF_forward_apont);
disp(polos_tf_apont);
polos_tf_estab = pole(TF_forward_estab);
disp(polos_tf_estab);

%Calculando os polos considerando a forma padrão de segunda ordem e os
%parâmetros ts(Acomodação), Máximo sobressinal(Mp) e 0 < p < 1(variação do
%tempo de acomodação)
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

%Polos que "satisfazem" os requisitos
polos_req_apont = [-sigma_apont+1i*wd_apont -sigma_apont-1i*wd_apont];
polos_req_estab = [-sigma_estab+1i*wd_estab -sigma_estab-1i*wd_estab];

disp(polos_req_apont);
disp(polos_req_estab);

%Considerando um controlador Compesador Lead-Lag Gc(s), ou seja, a função G(s) do ramo direto
%G(s) = K(s+ZL)/(s+PL) * A(s)*P(s)

%Calculando o 'a' utilizando a condição de fase do LGR

%thetaP1_estab = atan((wd_estab-imag(polos_tf_estab(1)))/(sigma_estab-real(polos_tf_estab(1))));
thetaP2_estab = atan((wd_estab-imag(polos_tf_estab(2)))/(sigma_estab-real(polos_tf_estab(2))));

PL = sigma_estab - wd_estab/tan(pi - thetaP2_estab);

%Calculando o 'K' utilizando a condição do módulo do LGR

%LP1_estab = sqrt((wd_estab-imag(polos_tf_estab(1)))^2 + (sigma_estab-real(polos_tf_estab(1)))^2);
LP2_estab = sqrt((wd_estab-imag(polos_tf_estab(2)))^2 + (sigma_estab-real(polos_tf_estab(2)))^2);

LPL_estab = sqrt((wd_estab-0)^2 + (sigma_estab-PL)^2);

K_estab = (LP2_estab*LPL_estab)/1.0;

disp(PL);
disp(K_estab);


function y = poles(Mp,p,ts)
    csi = sqrt(log(Mp)*log(Mp)/(pi^2 + log(Mp)^2));
    wn = (log(100) - log(100*p) - log(sqrt(1 - csi^2)))/(ts*csi);
    sigma = wn*csi;
    wd = wn*sqrt(1-csi^2);
    y = [wn,sigma,wd,csi];
end

