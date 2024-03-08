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

%Considerando um controlador PD, ou seja, a função G(s) do ramo direto
%G(s) = (Kp + N*kd) * (s+a)/(s+N) * A(s)*P(s)

%Calculando o 'a' utilizando a condição de fase do LGR

thetaP1_estab = atan((wd_estab-imag(polos_tf_estab(1)))/(sigma_estab-real(polos_tf_estab(1))));
thetaP2_estab = atan((wd_estab-imag(polos_tf_estab(2)))/(sigma_estab-real(polos_tf_estab(2))));

a_estab = abs(sigma_estab - wd_estab/tan(pi + thetaP1_estab + thetaP2_estab));


%Calculando o 'K' utilizando a condição do módulo do LGR

LP1_estab = sqrt((wd_estab-imag(polos_tf_estab(1)))^2 + (sigma_estab-real(polos_tf_estab(1)))^2);
LP2_estab = sqrt((wd_estab-imag(polos_tf_estab(2)))^2 + (sigma_estab-real(polos_tf_estab(2)))^2);
LP3_estab = sqrt((wd_estab-0)^2 + (sigma_estab-N)^2);

LZ1_estab = sqrt((wd_estab-0)^2 + (sigma_estab-a_estab)^2);

K_estab = (LP1_estab*LP2_estab*LP3_estab)/(LZ1_estab*K_tf_forward_estab);

Kd_estab = ((N-a_estab)*K_estab)/N^2;
Kp_estab = (a_estab*N*Kd_estab)/(N-a_estab);

thetaP1_apont = atan((wd_apont-imag(polos_tf_apont(1)))/(sigma_apont-real(polos_tf_apont(1))));
thetaP2_apont = atan((wd_apont-imag(polos_tf_apont(2)))/(sigma_apont-real(polos_tf_apont(2))));
thetaP3_apont = atan((wd_apont-imag(polos_tf_apont(3)))/(sigma_apont-real(polos_tf_apont(3))));
thetaP4_apont = atan((wd_apont-0)/(N-sigma_apont));
sum_theta_apont = (thetaP1_apont + thetaP2_apont+thetaP3_apont + thetaP4_apont);

a_apont = abs(sigma_apont - wd_apont/tan(-pi - pi - sum_theta_apont));
disp(a_apont);

LP1_apont = sqrt((wd_apont-imag(polos_tf_apont(1)))^2 + (sigma_apont-real(polos_tf_apont(1)))^2);
LP2_apont = sqrt((wd_apont-imag(polos_tf_apont(2)))^2 + (sigma_apont-real(polos_tf_apont(2)))^2);
LP3_apont = sqrt((wd_apont-imag(polos_tf_apont(3)))^2 + (sigma_apont-real(polos_tf_apont(3)))^2);
LP4_apont = sqrt((wd_apont-0)^2 + (sigma_apont-N)^2);

LZ1_apont = sqrt((wd_apont-0)^2 + (sigma_apont-a_apont)^2);

K_apont = (LP1_apont*LP2_apont*LP3_apont*LP4_apont)/(LZ1_apont*K_tf_forward_apont);

Kd_apont = ((N-a_apont)*K_apont)/N^2;
Kp_apont = (a_apont*N*Kd_apont)/(N-a_apont);

%Ajuste
%Kd_apont = Kd_apont*K_tf_forward_apont;
%Kp_apont = Kp_apont*K_tf_forward_apont;

function y = poles(Mp,p,ts)
    csi = sqrt(log(Mp)*log(Mp)/(pi^2 + log(Mp)^2));
    wn = (log(100) - log(100*p) - log(sqrt(1 - csi^2)))/(ts*csi);
    sigma = wn*csi;
    wd = wn*sqrt(1-csi^2);
    y = [wn,sigma,wd,csi];
end