clc;clear;format long g
%==========================================================================
%示例程序：
%1）如何将轨道根数转化为位置速度
%2）如何在J2摄动模型下进行轨道递推
%==========================================================================
%==========================================================================

global mu;
global J2;
global Re;

mu = 398600;  %地球引力常数 km2/s2
J2 = 0.0010826; %J2摄动常数
Re = 6.378137e3;       %地球半径 km 

%积分精度  设计结果检验时需要提供的积分精度为 1e-13 
RelTol = 1.e-13 ; AbsTol = 1.e-13; 

%积分时间 单位为秒
dt = 86400;
%轨道递推模型
MODEL = 'modelJ2';
TSPAN = [0 dt] ; 
%半长轴[km]     偏心率   轨道倾角[deg]  升交点经度[deg]  近地点幅角[deg]  真近点角[deg] 
ele_i = [6778.137    0.001    42.8   10   90   20];

%轨道根数转换为位置速度
rv_i = E2C(ele_i,mu);
tic
%数值积分
OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);
[t,rv] = ode113(MODEL,TSPAN,rv_i,OPTIONS) ;
toc
%输出最后一行的结果
rv_f = rv(end,:);

fprintf('初始时刻半长轴           =%e  [km]\n', ele_i(1))
fprintf('初始时刻偏心率           =%e  [-]\n', ele_i(2))
fprintf('初始时刻轨道倾角         =%e  [deg]\n', ele_i(3))
fprintf('初始时刻升交点赤经       =%e  [deg]\n', ele_i(4))
fprintf('初始时刻近地点幅角       =%e  [deg]\n', ele_i(5))
fprintf('初始时刻真近点角         =%e  [deg]\n\n', ele_i(6))

fprintf('轨道转移时间             =%e  [sec]\n', dt)
fprintf('初始时刻位置矢量         =[%e   %e  %e]  [km]\n', rv_i(1:3))
fprintf('初始时刻速度矢量         =[%e   %e  %e]  [km/s]\n', rv_i(4:6))
fprintf('初始时刻位置矢量         =[%e   %e  %e]  [km]\n', rv_f(1:3))
fprintf('初始时刻速度矢量         =[%e   %e  %e]  [km/s]\n', rv_f(4:6))




