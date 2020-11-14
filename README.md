# EKF-for-vehicel-state-estimation
6 states estimation of the vehicle fusing IMU and GPS
This project is about the state estimation of the vehicle with EKF under the co-simulation of the Simulink and Carsim.
The detail codes are included in the simulink file .mdl, and the configurations about the parameters of Carsim are included 
in the file of .par. You should import this file befor simulation.

The code of EKF as following:
function [X_est,Z_Cur] = fcn(U_Cur,X_Cur)
persistent X_est_1 Pk 
if isempty(X_est_1)
    X_est_1 = X_Cur;
    %% 误差协方差矩阵
Pk = 0.1*eye(6);
    %% 误差协方差矩阵

end

%% 常数设置
Ts = 1e-3;%跟仿真步长设置有关
Iz = 1343;%转动惯量
%% 系统误差向量
q = [2e-1 2e-1 1e-6 5e-1 1e-5 1e-5 ];
% q = zeros(1,6);
%% 测量误差向量
r = [0.2 0.2 0.00349 0.00349];
%% 输入
u1 = U_Cur(1);
u2 = U_Cur(2);
u3 = U_Cur(3)/Iz;
%% 测量加噪
Z_C1= X_Cur(1)+r(1)*(1-2*rand);
Z_C2= X_Cur(2)+r(2)*(1-2*rand);
Z_C3= X_Cur(3)+r(3)*(1-2*rand);
Z_C4= X_Cur(4)+r(4)*(1-2*rand);

Z_Cur =[Z_C1 Z_C2 Z_C3 Z_C4]';
%% 上一步最有估计状态
x1 = X_est_1(1);
x2 = X_est_1(2);
x3 = X_est_1(3);
x4 = X_est_1(4);
x5 = X_est_1(5);
x6 = X_est_1(6);

%% 噪声矩阵
Qk = diag(q)^2;
Rk = diag(r)^2;
%% 当前测量
% Z_Cur = Z_Ce_Curr;
%% 状态转移矩阵 
Ad=[1 0 -Ts*x5*sin(x3)-Ts*x6*cos(x3) 0 Ts*cos(x3) -Ts*sin(x3);
    0 1 Ts*x5*cos(x3)-Ts*x6*sin(x3) 0 Ts*sin(x3) Ts*cos(x3);
    0 0 1 Ts 0 0;
    0 0 0 1 0 0;
    0 0 0 Ts*x6 1 Ts*x4;
    0 0 0 -Ts*x5 -Ts*x4 1];
%% 观测矩阵
Cd=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0];
%% 模型估计值
x1_next = Ts*(x5*cos(x3)-x6*sin(x3))+x1;
x2_next = Ts*(x5*sin(x3)+x6*cos(x3))+x2;
x3_next = Ts*x4 + x3;
x4_next = Ts*u3 + x4;
x5_next = Ts*u1 +Ts*x4*x6 + x5;
x6_next = Ts*u2 - Ts*x4*x5 + x6;
x_forecast = [x1_next x2_next x3_next x4_next x5_next x6_next]';

%% 测量预测
z1 = x1_next;
z2 = x2_next;
z3 = x3_next;
z4 = x4_next;
Z_yuce = [z1 z2 z3 z4]';
%% 误差协方差矩阵预更新
Pkk_1 = Ad*Pk*Ad'+Ts^2*Qk;
%% 卡尔曼增益计算
%  Kk = Pkk_1*Cd'*(Cd*Pkk_1*Cd'+Rk)^-1;       %计算增益
 Kk = Pkk_1*Cd'*inv(Cd*Pkk_1*Cd'+Rk);       %计算增益
%% 最优估计
 x_hat = x_forecast+Kk*(Z_Cur-Z_yuce);      %校正
 %% 更新协方差矩阵
 Pk = (eye(6)-Kk*Cd)*Pkk_1;
 X_est = x_hat;
 X_est_1 = X_est;

Any problem please email to lubingev@sina.com.
