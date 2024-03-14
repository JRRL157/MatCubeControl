clc;
clear;
close all;
load('parametros4.mat');

% Equações de espaço de estados do "ramo direto"(Atuador + Dinâmica + Cinemática**)
A_estab = [-Ra/La -Kb/La 0;Ki/Jm -Bm/Jm 0;(Jr*Ki)/(Jb(1,1)*Jm) -(Jr*Bm)/(Jb(1,1)*Jm) 0];
B_estab = [1/La; 0; 0];

A_estab(3,1) = A_estab(3,1) + 1e-6;

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

%SIMULAÇÕES

simulation_time = 100;

% Simulink

out_sim = sim('controle_atitude_lqr',simulation_time);

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
        
        if mod(int32(time),5)==0
            set(q_ph1,'XData',time_arr(1:i),'YData',scalar_arr(1:i));
            set(q_ph2,'XData',time_arr(1:i),'YData',Q1_arr(1:i));
            set(q_ph3,'XData',time_arr(1:i),'YData',Q2_arr(1:i));
            set(q_ph4,'XData',time_arr(1:i),'YData',Q3_arr(1:i));
            
            set(vel_ph1,'XData',time_arr(1:i),'YData',omegaX_arr(1:i));
            set(vel_ph2,'XData',time_arr(1:i),'YData',omegaY_arr(1:i));
            set(vel_ph3,'XData',time_arr(1:i),'YData',omegaZ_arr(1:i));
            drawnow;
        end
        disp(time)
        %java.lang.Thread.sleep(0,time_step*1e4) % (ms, ns)
    end

    % plot(time_arr,scalar_arr);
end


