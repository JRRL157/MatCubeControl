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

%SIMULAÇÕES

simulation_time = 100;

% Simulink

out_sim = sim('controle_atitude',simulation_time);

scalar_arr = out_sim.simout;
Q1_arr = out_sim.simout1;
Q2_arr = out_sim.simout2;
Q3_arr = out_sim.simout3;
r1_theta_arr = out_sim.simout4;
r2_theta_arr = out_sim.simout5;
r3_theta_arr = out_sim.simout6;
omegaX_arr = out_sim.simout7;
omegaY_arr = out_sim.simout8;
omegaZ_arr = out_sim.simout9;
time_arr = out_sim.tout;

% CoppeliaSim

vrep = remApi('remoteApi');
vrep.simxFinish(-1);

clientId = vrep.simxStart('127.0.0.1',19999,true,true,5000,5);

handle_ret = 0;
num_ite = 0;

if clientId < 0
    disp("Connection FAILED!!!");
    vrep.delete;
    return;
else
    fprintf("Connection %d to remote API server\n",clientId);
    
    [ret, rotor_x] = vrep.simxGetObjectHandle(clientId,'motor_x',vrep.simx_opmode_oneshot_wait);
    [ret, rotor_y] = vrep.simxGetObjectHandle(clientId,'motor_y',vrep.simx_opmode_oneshot_wait);
    [ret, rotor_z] = vrep.simxGetObjectHandle(clientId,'motor_z',vrep.simx_opmode_oneshot_wait);
    [ret, cubesat] = vrep.simxGetObjectHandle(clientId,'cubesat',vrep.simx_opmode_oneshot_wait);    
    time = 0;
    
    %Plots in Real time
    fig1 = figure(1);
    ax1 = axes('NextPlot','add');       
    legend(ax1,{'escalar','v1','v2','v3'});
    q_ph1 = plot(ax1,nan,nan);    
    q_ph2 = plot(ax1,nan,nan);            
    q_ph3 = plot(ax1,nan,nan);            
    q_ph4 = plot(ax1,nan,nan);
    
    
    fig2 = figure(2);
    ax2 = axes('NextPlot','add');
    vel_ph1 = plot(ax2,nan,nan);
    vel_ph2 = plot(ax2,nan,nan);
    vel_ph3 = plot(ax2,nan,nan);

    for i=2:length(scalar_arr)
        vrep.simxSetObjectQuaternion(clientId,cubesat,-1,[scalar_arr(i-1),Q1_arr(i-1),Q2_arr(i-1),Q3_arr(i-1)],vrep.simx_opmode_streaming);
        vrep.simxSetJointTargetPosition(clientId,rotor_x,r1_theta_arr(i-1),vrep.simx_opmode_streaming);
        vrep.simxSetJointTargetPosition(clientId,rotor_y,r2_theta_arr(i-1),vrep.simx_opmode_streaming);
        vrep.simxSetJointTargetPosition(clientId,rotor_z,r3_theta_arr(i-1),vrep.simx_opmode_streaming);
        
        time_step = time_arr(i)-time_arr(i-1);  % segundos
        time = time + time_step;
        
        if mod(num_ite,10000)==0
            if i <= 100000
                start_idx = 1;
            else
                start_idx = i - 100000+5;
            end

            set(q_ph1,'XData',time_arr(start_idx:i),'YData',scalar_arr(start_idx:i));
            set(q_ph2,'XData',time_arr(start_idx:i),'YData',Q1_arr(start_idx:i));
            set(q_ph3,'XData',time_arr(start_idx:i),'YData',Q2_arr(start_idx:i));
            set(q_ph4,'XData',time_arr(start_idx:i),'YData',Q3_arr(start_idx:i));
            
            set(vel_ph1,'XData',time_arr(start_idx:i),'YData',omegaX_arr(start_idx:i));
            set(vel_ph2,'XData',time_arr(start_idx:i),'YData',omegaY_arr(start_idx:i));
            set(vel_ph3,'XData',time_arr(start_idx:i),'YData',omegaZ_arr(start_idx:i));
            drawnow;
        end
        disp(time)
        num_ite = num_ite + 1;
        %java.lang.Thread.sleep(0,time_step*1e4) % (ms, ns)
    end

    set(q_ph1,'XData',time_arr(1:i),'YData',scalar_arr(1:i));
    set(q_ph2,'XData',time_arr(1:i),'YData',Q1_arr(1:i));
    set(q_ph3,'XData',time_arr(1:i),'YData',Q2_arr(1:i));
    set(q_ph4,'XData',time_arr(1:i),'YData',Q3_arr(1:i));
    
    set(vel_ph1,'XData',time_arr(1:i),'YData',omegaX_arr(1:i));
    set(vel_ph2,'XData',time_arr(1:i),'YData',omegaY_arr(1:i));
    set(vel_ph3,'XData',time_arr(1:i),'YData',omegaZ_arr(1:i));
    drawnow;
end


function y = poles(Mp,p,ts)
    csi = sqrt(log(Mp)*log(Mp)/(pi^2 + log(Mp)^2));
    wn = (log(100) - log(100*p) - log(sqrt(1 - csi^2)))/(ts*csi);
    sigma = wn*csi;
    wd = wn*sqrt(1-csi^2);
    y = [wn,sigma,wd,csi];
end