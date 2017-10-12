function GainComp = Gain(TaskNo, TN, RVN1, RVN2, RVN3, lamda, El, Az) 
%==========================================================================
%   曾豪   zenghaohz@163.com      2017-09-10
%==========================================================================
% 功能：求当前子任务中对GEO卫星的监测收益
% 输入： 
%      TaskNo         任务编号
%      TN             整秒时间序列，N维列向量，单位：s
%      RVN1           波束内对应于TN的监测卫星1的状态序列，N*6维矩阵，单位：km，km/s
%      RVN2           波束内对应于TN的监测卫星2的状态序列，N*6维矩阵，单位：km，km/s
%      RVN3           波束内对应于TN的监测卫星3的状态序列，N*6维矩阵，单位：km，km/s
%      lamda          GEO卫星定点经度，deg
%      El             波束指向俯仰角，deg
%      Az             波束指向方位角，deg
% 输出：
%      GainComp       当前子任务中对GEO卫星的监测收益，单位：km^2
%==========================================================================

RN = [RVN1(:,1:3)'; RVN2(:,1:3)'; RVN3(:,1:3)'];
GainN = zeros(length(TN),1);

% 参数
rGEO = 42164.2; % GEO轨道半径，单位：km
Thetag0 = 100.1218 * pi / 180; % 初始时刻(MJD = 58849)的格林尼治恒星时角，单位：rad
Dthetag = 360 / 0.99727 / 86400 * pi / 180; % 格林尼治恒星时角的变化率，单位：rad/sec

% 单位换算：deg -> rad
lamda = lamda * pi / 180;
El = El * pi / 180;
Az = Az * pi / 180;

for jj = 1:length(TN)
    
    % GEO卫星在J2000坐标系下的位置矢量
    Thetag = Thetag0 + Dthetag * TN(jj); %rad
    Mie = [ cos(-Thetag)  sin(-Thetag) 0;
           -sin(-Thetag)  cos(-Thetag) 0;
            0            0           1];
    rGEO_ECF = [rGEO * cos(lamda); rGEO * sin(lamda); 0];
    rGEO_J2000 = Mie * rGEO_ECF;

    % GEO卫星波束矢量在J2000坐标系下的表示
    g = [cos(El) * cos(Az); cos(El) * sin(Az); sin(El)]; % % 波束矢量（轨道坐标系）
    Mx90 = [1,  0,               0;
            0,  cos(90*pi/180),  sin(90*pi/180);
            0, -sin(90*pi/180),  cos(90*pi/180)];
    Mz90 = [ cos(-90*pi/180), sin(-90*pi/180), 0;
            -sin(-90*pi/180), cos(-90*pi/180), 0;
             0,               0,               1];
    MzTL = [ cos(- Thetag - lamda), sin(- Thetag - lamda), 0;
            -sin(- Thetag - lamda), cos(- Thetag - lamda), 0;
             0,                     0,                     1];
    M     = MzTL * Mz90 * Mx90;
    Beam = M * g; % 波束矢量（J2000坐标系）
    
    for kk = 1 : 3
        
        % 监测卫星相对于GEO的位置矢量
        Ri = RN(1+(kk-1)*3 : 3+(kk-1)*3,jj);
        Rgi = Ri - rGEO_J2000;

        % 检测监测卫星是否在波束内
        Angle1 = acos( dot(Rgi, Beam) / (norm(Rgi) * norm(Beam)));
%         if Angle1 * 180 / pi > 0.5 
%             fprintf('第 %d 次子任务中 第 %d 颗监测卫星不在波束圆锥内\n',TaskNo, kk);
%         end
        
        % 检测监测卫星与GEO卫星是否在地球同侧 
        Angle2 = acos( dot(Ri, rGEO_J2000) / (norm(Ri) * norm(rGEO_J2000)));
%         if Angle2 * 180 / pi  > 90
%             fprintf('第 %d 次子任务中 第 %d 颗监测卫星与GEO卫星不在同侧\n',TaskNo, kk);
%         end        
    end
    
    % 计算监测收益
    R12 = RN(1:3,jj) - RN(4:6,jj);
    R23 = RN(4:6,jj) - RN(7:9,jj);
    R13 = RN(1:3,jj) - RN(7:9,jj);
    L = (norm(R12) + norm(R13) + norm(R23))/2;
    S = sqrt( L * (L - norm(R12)) * ( L- norm(R13)) * (L - norm(R23)) );
    temp = cross(R12, R13);
    gk = abs(dot(Beam, temp/norm(temp))); % 波束与监测卫星平面夹角的正弦
    GainN(jj) = S * gk;     
end

% 当前子任务的监测总收益
GainComp = sum(GainN);

end

