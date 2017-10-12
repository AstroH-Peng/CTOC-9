%==========================================================================
%               第九届全国空间轨道设计竞赛甲组结果验算程序
%==========================================================================
%   祁瑞     qirui@bit.edu.cn      2017-09-10
%==========================================================================

%清空工作区
clc;clear

%==========================================================================
%精度不满足要求，根据情况判定为轻度违反约束和重度违反约束，判定原则由命题方决定
%将给2天的时间进行修改，修改后得分不能超过原方案得分，飞行序列须与原方案基本一致
%==========================================================================

%% =========================参数设置===============================
%初始化全局变量
global mu J2 Re ge Isp 
global rTol vTol mTol

mu   = 398600;    %地球引力常数 km2/s2
J2   = 0.0010826; %J2摄动常数
Re   = 6378.137;  %地球半径 km
ge   = 9.80665*1e-3;   %地球海平面的重力加速度， km/sec^2
Isp  = 300;       %推进器比冲，单位：s
SPstep = 1;       %检测星间距离的时间采样步长，单位：s

% 位置、速度、质量、收益验算精度
rTol = 1e-3; % 单位：km
vTol = 1e-6; % 单位：km/s
mTol = 1e-4; % 单位：kg
GainTol = 1e-3; % 单位：km^2

% 积分器设置误差和步长
RelTol = 1e-13 ; 
AbsTol = 1e-13;
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

GEOgain = zeros(18,1);

%% ==========================数据读入==============================
% 参赛者数据
data = load('result_test.txt');

% 官方GEO卫星经度、波束俯仰角、方位角
GEO = load('GEO.txt');


%% ========================初始状态检测============================
% 初始轨道要素和质量
sma0 = 6778.137; %半长轴，单位：km
ecc0 = 0; %偏心率
inc0 = 42.8*pi/180; %轨道倾角，单位：rad
m0 = 500; %监测卫星质量，单位：kg

% 检验三颗监测卫星的初始根数误差
fprintf('================================初始状态检测==================================== \n\n')
for kk = 1:3
    flagIni = 0;
    rv0(kk,:)  = data(min(find(data(:,4) == kk)),5:10);
    ele0(kk,:) = RV2Ele(rv0(kk,:),mu);
    dsma = ele0(kk,1) - sma0;
    decc = ele0(kk,2) - ecc0;
    dinc = ele0(kk,3) - inc0;
    dm   = data(min(find(data(:,4) == kk)),11) - m0;
    if abs(dsma) > 1e-5
        fprintf('第%d颗监测卫星违反初始半长轴约束：名义半长轴%f，实际半长轴%f\n\n',kk, sma0, ele0(kk,1))
        flagIni = flagIni + 1;
    end
    if abs(decc) > 1e-7
        fprintf('第%d颗监测卫星违反初始偏心率约束：名义偏心率%f，实际偏心率%f\n\n',kk, ecc0, ele0(kk,2))
        flagIni = flagIni + 1;
    end
    if abs(dinc) > 1e-7
        fprintf('第%d颗监测卫星违反初始轨道倾角约束：名义倾角%f，实际倾角%f\n\n',kk, inc0, ele0(kk,3))
        flagIni = flagIni + 1;
    end
    if abs(dm) > 1e-10
        fprintf('第%d颗监测卫星违反初始质量约束：名义质量%f，实际质量%f\n\n',kk, m0, data(min(find(data(:,4) == kk)),11))
        flagIni = flagIni + 1;
    end
    if flagIni == 0
        fprintf('第%d颗监测卫星初始状态验证成功 \n\n', kk)
    else
        fprintf('第%d颗监测卫星初始状态验证失败 \n\n', kk)
    end
end
fprintf('=============================================================================== \n\n')

%% =============轨道积分，精度检测，提取监测信息，计算收益=============
% 初始状态
T10 = 0;
T20 = 0;
T30 = 0;
T1End = -1;
T2End = -1;
T3End = -1;
RV10 = rv0(1,:);
RV20 = rv0(2,:);
RV30 = rv0(3,:);
m10 = m0;
m20 = m0;
m30 = m0;

for ii = 1 : data(end,2) % data(end,2)是子任务个数
    
    % 数据读取
    SubTask =  data(find(data(:,2) == ii),:); % 读取子任务ii的所有数据，赋给SubTask
    Sat1 = SubTask(find(SubTask(:,4) == 1),:); % 读取子任务ii中各监测卫星的所有数据，分别赋给Sat1、Sat2、Sat3
    Sat2 = SubTask(find(SubTask(:,4) == 2),:);
    Sat3 = SubTask(find(SubTask(:,4) == 3),:);
    
    if Sat1(find(Sat1(:,15)==1),3) == T1End
        disp('监测卫星1在上一个子任务出波束时间不可与当前子任务进波束时间相同')
    end
    if Sat2(find(Sat2(:,15)==1),3) == T2End
        disp('监测卫星2在上一个子任务出波束时间不可与当前子任务进波束时间相同')
    end
    if Sat3(find(Sat3(:,15)==1),3) == T3End
        disp('监测卫星3在上一个子任务出波束时间不可与当前子任务进波束时间相同')
    end
    
    %Orbit函数实现分段轨道积分，检测位置、速度、质量精度和轨道高度，提取轨道和监测状态点信息
    [TH1, RVH1, Tobv10, Tobv1f, RVobv10] = Orbit(T10, RV10, m10, Sat1, options, ii, 1);
    if any((sqrt(RVH1(:,1).^2 + RVH1(:,2).^2 + RVH1(:,3).^2) - Re) <300) == 1
        fprintf('第 %d 个子任务中 卫星 %d 轨道高度小于 300 km\n', ii, 1);
    end
    
    [TH2, RVH2, Tobv20, Tobv2f, RVobv20] = Orbit(T20, RV20, m20, Sat2, options, ii, 2);
    if any((sqrt(RVH2(:,1).^2 + RVH2(:,2).^2 + RVH2(:,3).^2) - Re) <300) == 1
        fprintf('第 %d 个子任务中 卫星 %d 轨道高度小于 300 km\n', ii, 2);
    end
    
    [TH3, RVH3, Tobv30, Tobv3f, RVobv30] = Orbit(T30, RV30, m30, Sat3, options, ii, 3);
    if any((sqrt(RVH3(:,1).^2 + RVH3(:,2).^2 + RVH3(:,3).^2) - Re) <300) == 1
        fprintf('第 %d 个子任务中 卫星 %d 轨道高度小于 300 km\n', ii, 1);
    end
    
    % 检验各监测卫星进出波束时间的一致性
    if Tobv10 ~= Tobv20 || Tobv10 ~= Tobv30
        fprintf('第 %d 个子任务中三颗监测卫星进波束的时间不一致\n', ii);
    end
    if Tobv1f ~= Tobv2f || Tobv1f ~= Tobv3f
        fprintf('第 %d 个子任务中三颗监测卫星出波束的时间不一致\n', ii);
    end
    
    % 检测三颗监测卫星间的相对距离是否满足约束
    Tsp = TH1(1) : SPstep : TH1(end);
    RVsp1 = spline(TH1', RVH1', Tsp);
    RVsp2 = spline(TH2', RVH2', Tsp);
    RVsp3 = spline(TH3', RVH3', Tsp);
    for gg = 1 : length(Tsp)
        if norm(RVsp1(1:3, gg) - RVsp2(1:3, gg)) < 0.1
            fprintf('第 %d 个子任务中第 1 和 第 2 颗监测卫星之间的距离在 t = %f 时小于100m\n', ii, Tsp(gg));
        elseif norm(RVsp2(1:3, gg) - RVsp3(1:3, gg)) < 0.1
            fprintf('第 %d 个子任务中第 2 和 第 3 颗监测卫星之间的距离在 t = %f 时小于100m\n', ii, Tsp(gg));
        elseif norm(RVsp1(1:3, gg) - RVsp3(1:3, gg)) < 0.1
            fprintf('第 %d 个子任务中第 1 和 第 3 颗监测卫星之间的距离在 t = %f 时小于100m\n', ii, Tsp(gg));
        end
    end
 
    % 提取波束内各监测卫星在整秒时的状态     
    tspan = [Tobv10, ceil(Tobv10) : 1 : Tobv1f];
    flag = 0;        
    if tspan(1) == tspan(2)
        tspan = tspan(2:end);
        flag = 1;
    end
        
    [tN1,yN1]  = ode113('modelJ2', tspan, RVobv10, options) ;  
    [~,yN2]    = ode113('modelJ2', tspan, RVobv20, options) ;  
    [~,yN3]    = ode113('modelJ2', tspan, RVobv30, options) ;  
    
    if size(tspan,2) == 2
        ttN1 = [tN1(1); tN1(end)];
        yyN1 = [yN1(1,:); yN1(end,:)];
        yyN2 = [yN2(1,:); yN2(end,:)];
        yyN3 = [yN3(1,:); yN3(end,:)]; 
        
    elseif size(tspan,2) == 1
        ttN1 = tspan;
        yyN1 = RVobv10;
        yyN2 = RVobv20;
        yyN3 = RVobv30;         
        
    elseif flag == 0
        ttN1 = tN1(2:end);
        yyN1 = yN1(2:end,:);
        yyN2 = yN2(2:end,:);
        yyN3 = yN3(2:end,:);
    else
        ttN1 = tN1;
        yyN1 = yN1;
        yyN2 = yN2;
        yyN3 = yN3;        
    end
    
    % 计算收益，检测各监测卫星均在波束内，且在地球同侧
    GEONo  = Sat1(end,16);
    lamda  = GEO(GEONo,2);
    El     = GEO(GEONo,3);
    Az     = GEO(GEONo,4);    
    GainComp = Gain(ii, ttN1, yyN1, yyN2, yyN3, lamda, El, Az); % 计算出的当前子任务获得的GEO监测收益  
    GainProv = (SubTask(end,17) - SubTask(1,17)); % 竞赛者提供的当前子任务获得的GEO监测收益  
    if abs(GainComp - GainProv) > GainTol
        fprintf('第 %d 个子任务的监测收益不正确,收益绝对误差为 %8f (km^2)\n',ii,abs(GainComp - GainProv));
    end
    GEOgain(GEONo) = GEOgain(GEONo) + GainProv;
    
    % 当前子任务的末状态赋值为下一个子任务的初状态

    T1End = Sat1(find(Sat1(:,15)==2),3);
    T2End = Sat2(find(Sat2(:,15)==2),3);
    T3End = Sat3(find(Sat3(:,15)==2),3);
    T10 = TH1(end);
    T20 = TH2(end);
    T30 = TH3(end);
    RV10 = RVH1(end,:);
    RV20 = RVH2(end,:);
    RV30 = RVH3(end,:);
    m10 = Sat1(end,11);
    m20 = Sat2(end,11);
    m30 = Sat3(end,11);    


    fprintf('第 %d / %d 个子任务验证结束 \n\n', ii, data(end,2)) % data(end,2)是子任务个数
end


%% =========================结果输出===============================
if min(GEOgain) >= 1e7
    Obj1 = Tobv1f;
    mf = [Sat1(end,11);Sat2(end,11);Sat3(end,11)];
    if min(m0-mf) > 200
        error('燃料消耗超过 200 kg');
    end
    Obj2 = sum(m0 - mf);
    disp('总任务已完成！')
    
    % 输出设计指标
    disp('第一指标 (day)')
    disp(Obj1/3600/24) % 输出第一指标，单位：day
    disp('第二指标 (kg)')
    disp(Obj2)
    disp('各监测卫星燃料消耗 (kg)')
    disp(m0 - mf)
        
    % 输出各GEO卫星的收益
    disp('GEO卫星收益')
    disp(GEOgain)

else
    disp('总任务未完成！')
end