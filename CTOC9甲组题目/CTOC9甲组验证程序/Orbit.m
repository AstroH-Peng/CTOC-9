function [TH, RVH, Tobv0, Tobvf, RVobv0] = Orbit(T0, RV0, m0, Sat, options, TaskNo, SatNo)
%==========================================================================
%   杨科莹    yangkeying@bit.edu.cn      2017-09-10
%==========================================================================
% 功能：实现分段轨道积分，检测位置、速度、质量精度，提取轨道和监测状态点信息
% 输入： 
%      T0           子任务初始时刻，单位：s
%      RV0          监测卫星的初始状态，6维列向量，单位：km, km/s
%      m0           监测卫星的初始质量，单位：kg
%      Sat          子任务中该监测卫星的数据矩阵，N*17维矩阵
%      options      积分器设置
%      TaskNo       当前子任务序号
%      SatNo        该监测卫星序号
% 输出：
%      TH           当前子任务中该监测卫星的时间采样点，单位：s
%      RVH          当前子任务中该监测卫星状态的时间历程，单位：km, km/s
%      Tobv0        该监测卫星进入监测的时刻，单位：s
%      Tobvf        该监测卫星结束监测的时刻，单位：s
%      RVobv0       该监测卫星进入监测时的状态，6维列向量，单位：km, km/s
%==========================================================================

global Re mu Isp ge

TH = [];
RVH = [];
Count1 = 0; % 用于记录状态1出现的次数
Count2 = 0; % 用于记录状态2出现的次数
Count3 = 0; % 用于记录状态3出现的次数

% 检测时间单调性
Ttemp = Sat(:,3);
if any(diff(Ttemp)<0)
    fprintf('第 %d 个子任务的时间不是单调递增\n',TaskNo);
end

% 逐行积分
for jj = 1 : size(Sat,1)    
    tspan = [T0,Sat(jj,3)];
    if (tspan(2) - tspan(1) == 0)  && (Sat(jj,15) == 1)
        Tobv0  = Sat(jj,3);
        RVobv0 = Sat(jj,5:10);
        continue;
    elseif (tspan(2) - tspan(1) == 0)  && (Sat(jj,15) == 3)
        Tobv0  = Sat(jj,3);
        RVobv0 = Sat(jj,5:10);
        deltaV = norm(Sat(jj,12:14));
        deltam = m0 * (1- exp( - deltaV / (Isp * ge)));
        mf = m0 - deltam;
        testRVM(TaskNo, SatNo, jj, RV0(1:3), Sat(jj,5:7), RV0(4:6), Sat(jj,8:10) - Sat(jj,12:14), mf, Sat(jj,11));
        Count3 = Count3 + 1; 
        
        T0 = Sat(jj,3);
        RV0 = Sat(jj,5:10);
        m0 = Sat(jj,11);
        
        continue;
    elseif (tspan(2) - tspan(1) == 0)
        Tobvf = Sat(jj,3);
        continue;
    end
    [tt,yy] = ode113('modelJ2', tspan, RV0, options);
    
    % ==================  eventID = 4  ====================================
    if Sat(jj,15) == 4
        ;
    % ==================  eventID = 3  ====================================
    elseif Sat(jj,15) == 3
        deltaV = norm(Sat(jj,12:14));
        deltam = m0 * (1- exp( - deltaV / (Isp * ge)));
        mf = m0 - deltam;
        testRVM(TaskNo, SatNo, jj, yy(end,1:3), Sat(jj,5:7), yy(end,4:6), Sat(jj,8:10) - Sat(jj,12:14), mf, Sat(jj,11));
        Count3 = Count3 + 1;  
    % ==================  eventID = 1  ====================================  
    elseif Sat(jj,15) == 1
        testRVM(TaskNo, SatNo, jj, yy(end,1:3), Sat(jj,5:7), yy(end,4:6), Sat(jj,8:10), m0, Sat(jj,11));
        Tobv0 = Sat(jj,3);
        RVobv0 = Sat(jj,5:10);
        Count1 = Count1 + 1;
        
    % ==================  eventID = 2  ====================================
    elseif Sat(jj,15) == 2        
        testRVM(TaskNo, SatNo, jj, yy(end,1:3), Sat(jj,5:7), yy(end,4:6), Sat(jj,8:10), m0, Sat(jj,11));
        Tobvf = Sat(jj,3);
        Count2 = Count2 + 1;
    else
        fprintf('event的取值只能是 1, 2, 3, 4.\n');
    end
    
    % 将当前轨道段的末状态赋值为下一轨道段的初状态
    T0 = Sat(jj,3);
    RV0 = Sat(jj,5:10);
    m0 = Sat(jj,11);
    
    % 记录整段子任务轨道的时刻和状态时间历程
    TH = cat(1, TH, tt(1:end-1));
    RVH = cat(1, RVH, yy(1:end-1,1:6));
end
if Count1 > 1 || Sat(end-1,15) ~= 1
    fprintf('第 %d 个子任务中 卫星 %d 中状态1出现的次数或位置有误\n',TaskNo, SatNo);
end
if Count2 > 1 || Sat(end,15) ~= 2
    fprintf('第 %d 个子任务中 卫星 %d 中状态2出现的次数或位置有误\n',TaskNo, SatNo);
end
if Count3 > 5
    fprintf('第 %d 个子任务中 卫星 %d 的机动次数超过了 5 次\n',TaskNo, SatNo);
end


function testRVM(TaskNo, SatNo, ColNo, rf, rt, vf, vt, mf, mt)
% 用于检测位置、速度、质量精度

global rTol vTol mTol

if abs(norm(rf)-norm(rt)) > rTol
    fprintf('在第 %d 个子任务中，第 %d 颗监测卫星在第 %d 行的位置精度不符合要求,位置误差为 %8f (km)\n',TaskNo, SatNo, ColNo,abs(norm(rf)-norm(rt)))
end
if abs(norm(vf)-norm(vt)) > vTol
    fprintf('在第 %d 个子任务中，第 %d 颗监测卫星在第 %d 行的速度精度不符合要求,速度误差为 %8f (km/s)\n',TaskNo, SatNo, ColNo,abs(norm(vf)-norm(vt)))
end
if abs(norm(mf)-norm(mt)) > mTol
    fprintf('在第 %d 个子任务中，第 %d 颗监测卫星在第 %d 行的质量精度不符合要求,质量误差为 %8f (km)\n',TaskNo, SatNo, ColNo,abs(norm(mf)-norm(mt)))
end