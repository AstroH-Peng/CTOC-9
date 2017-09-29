% =========================================================================
% 第9届全国空间轨道设计竞赛乙组题目验证程序
% 对'constellation.txt'和'transfer.txt'进行检查
%
% 参考:
% 李明涛, 第8届全国空间轨道设计竞赛甲组题目验证程序
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 1.1.1.20170901
% -------------------------------------------------------------------------
% Verification program of problem B of CTOC9
% check 'constellation.txt' and 'transfer.txt'
%
% Reference:
% MingTao LI, the verification program of problem A of CTOC8
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% 1.1.1.Sep/01/2017
% -------------------------------------------------------------------------
function check
clc; clear; close all;

% 初始化参数 (initialize parameters)
aux = aux_;

% 选择提示语言(select language)
aux.langFlag = 1; % 中文
% aux.langFlag = 2; % English

if aux.langFlag == 1
    disp('==============================================================');
    disp('(1) 正在检查 constellation.txt')
else
    disp('==============================================================');
    disp('(1) checking constellation.txt')
end

% 读取 city.txt (Read city.txt)
[id , lon , lat , weight] = textread('city.txt' , '%s%f%f%f' , 'headerlines' , 1);
city.id = id;
city.lonlat = [lon , lat];
city.weight = weight;
city.ecf = [aux.req * cosd(lat) .* cosd(lon), ...
    aux.req * cosd(lat) .* sind(lon), ...
    aux.req * sind(lat)];
aux.city = city;

% 读取 carry.txt (Read carry.txt)
[id , MJD2000 , smaAray , eccAray , inclAray , raanAray , argperAray , manomAray] = ...
    textread('carry.txt' , '%s%f%f%f%f%f%f%f' , 'headerlines' , 1);
tempMtx = [smaAray , eccAray , inclAray , raanAray , argperAray , manomAray];
for i = 1 : size(id , 1)
    carry.id(i , :) = str2num(cell2mat(id(i)));
end
carry.MJD2000 = MJD2000;
carry.oe = tempMtx;
aux.carry = carry;

% 读取 constellation.txt (Read constellation.txt)
[id , MJD2000 , smaAray , eccAray , inclAray , raanAray , argperAray , manomAray] = ...
    textread('constellation.txt' , '%s%f%f%f%f%f%f%f' , 'headerlines' , 0);
tempMtx = [smaAray , eccAray , inclAray , raanAray , argperAray , manomAray];
for i = 1 : size(id , 1)
    cons.id(i , :) = str2num(cell2mat(id(i)));
end
cons.MJD2000 = MJD2000;
cons.oe = tempMtx;
aux.cons = cons;

% 检查：星座编号 (check the indexes of the constellation)
if (cons.id(1) ~= 1) || any(diff(cons.id) - 1)
    if aux.langFlag == 1
        disp('错误：星座第一颗卫星编号不为1或者编号不连续！')
    else
        disp('Warning: the first index is not 1 or the constellation indexes are not continuous!')
    end
    aux.errorFlag = 1;
end

% 检查：星座历元 (check the epoch of the constellation)
if any(MJD2000 - aux.MJD2000_tf)
    if aux.langFlag == 1
        disp('错误：星座开始服务时间错误！')
    else
        disp('Warning: the epoch that constellation begins service is wrong!')
    end
    aux.errorFlag = 1;
end

% 检查：星座轨道高度约束 (check the perigee altitude of the constellation)
if any(smaAray .* (1 - eccAray) - aux.req < aux.altiMin)
    if aux.langFlag == 1
        disp('错误：星座近地点高度低于500km！')
    else
        disp('Warning: the perigee altitude of the constellation is lower than 500km!')
    end
    aux.errorFlag = 1;
end

% 检查：星座椭圆轨道倾角约束 (check the inclination of the constellation)
if any( min( abs( [inclAray(eccAray ~= 0) - acos(sqrt(1 / 5)) , inclAray(eccAray ~= 0) - pi + acos(sqrt(1 / 5))] )' ) > 1e-6)
    if aux.langFlag == 1
        disp('错误：星座的椭圆轨道没有采用临界倾角！')
    else
        disp('Waring: the critical inclination is not used in elliptic orbits!')
    end
    aux.errorFlag = 1;
end

% 计算Obj1 (evaluate Obj1)
Obj1 = getObj1_(aux);

% 显示检查结果 (display the verification result of constellation.txt)
if aux.errorFlag == 0
    if aux.langFlag == 1
        disp('constellation.txt 通过验证！');
    else
        disp('constellation.txt passed the verification!');
    end
else
    if aux.langFlag == 1
        disp('constellation.txt 未能通过验证，请检查相关错误！')
    else
        disp('constellation.txt failed to pass the verification! please check errors!')
    end
end

if aux.langFlag == 1
    disp('==============================================================');
    disp('(2) 正在检查 transfer.txt')
else
    disp('==============================================================');
    disp('(2) checking transfer.txt')
end

visitAray = ones(size(cons.id)); % 统计入轨卫星，每入轨一颗卫星相应位置设为0 (count the number of injected satellites)
carryAray = zeros(size(carry.id)); % 统计使用的搭载任务 (count the number of piggyback missions)
NTask = 0; % 统计任务数 (count the number of tasks)
NsubTask = 0; % 统计子任务数 (count the number of subtasks)
Nnavi = 0; % 统计导航卫星数 (count the number of navigation satellites)
Nlaunch = 0; % 统计火箭数 (count the number of launching rockets)
Ncarry = 0; % 统计搭载数 (count the number of piggyback missions)

fid = fopen(['transfer.txt'] , 'r'); % 读取 'transfer.txt' (read 'transfer.txt')
tline = fgets(fid); % 读取第一行 (read the first line)
while tline ~= -1 % 读取'transfer.txt'直到最后一行 (read 'transfer.txt' until the end line)
    
    % 读取tline (read tline)
    lineTemp = textscan(tline ,'%s %f');
    Flag = lineTemp{1};
    
    % 如果检测到'Task' (if 'task' is detected)
    if strcmp(Flag , 'Task')
        
        % 主任务数 + 1 (number of task + 1)
        NTask = NTask + 1;
        
        % 子任务数清零 (clear subtask)
        NsubTask = 0;
        
        % 读取任务类型 (read tline)
        tline = fgets(fid);
        lineTemp = textscan(tline , '%s %s');
        
        % 识别任务类型 (identify the type of task)
        if strcmp(lineTemp{2} , 'Launch')
            Nlaunch = Nlaunch + 1; % 发射火箭数 + 1 (number of new launchings + 1)
            typeFlag = 1;
        else
            Ncarry = Ncarry + 1; % 搭载数 + 1 (number of piggyback missions + 1)
            typeFlag = 2;
        end
        
        % 读取停泊轨道状态 (read state of parking orbit)
        tline = fgets(fid);
        tline = fgets(fid);
        lineTemp = cell2mat(textscan(tline , '%f'))';
        
        % 构造停泊轨道结构体 (construct a structure of parking orbit)
        park.MJD2000 = lineTemp(1); % 停泊轨道历元 (epoch of parking orbit)
        park.oe = lineTemp(2:7); % 停泊轨道轨道根数 (orbital elements of parking orbit)
        [eciTemp , ~] = oe2eciecf_(park.oe , park.MJD2000 , park.MJD2000 , aux); % 轨道根数转位置速度 (from orbital elements to position and velocity)
        park.sv = eciTemp'; % 停泊轨道位置速度 (save position and velocity of parking orbit)
        
        % 检查：停泊轨道开始时间是否在建设时间约束内 (check: whether the epoch of parking orbit is located in the construction phase)
        if park.MJD2000 < aux.MJD2000_ti || park.MJD2000  > aux.MJD2000_tf
            if aux.langFlag == 1
                disp('错误：停泊轨道开始时间不在星座建设时间约束内！');
            else
                disp('Warning: the epoch of parking orbit is out of the construction phase!');
            end
            aux.errorFlag = 1;
        end
        
        % 检查：停泊轨道每颗卫星质量是否等于100kg (check: whether the mass of each satellite at parking orbit is equal to 100kg)
        if lineTemp(8) - aux.m0 ~= 0
            if aux.langFlag == 1
                disp('错误：停泊轨道单颗卫星质量不等于100kg！');
            else
                disp('Warning：the initiall mass of each satellite at parking orbit is not equal to 100kg!');
            end
            aux.errorFlag = 1;
        end
        
        % 检查停泊轨道状态 (check the state of parking orbit)
        switch typeFlag
            case 1 % 发射新火箭 (for new launchings)
                % 检查：停泊轨道高度是否为900km (check: whether the altitude of parking orbit is equal to 900km)
                if abs(park.oe(1) - (aux.altiPark + aux.req)) > aux.diffRMax
                    if aux.langFlag == 1
                        disp('错误：停泊轨道高度违反约束！')
                    else
                        disp('Warning: the altitude of parking orbit is not equal to 900km!')
                    end
                    aux.errorFlag = 1;
                end
                % 检查：停泊轨道是否为圆轨道 (check: whether the parking orbit is circular)
                if park.oe(2) > 1e-6
                    if aux.langFlag == 1
                        disp('错误：停泊轨道非圆轨道！')
                    else
                        disp('Warning: the parking orbit is not circular!')
                    end
                    aux.errorFlag = 1;
                end
                
            case 2 % 搭载 (for piggyback missions)
                % 计算所有搭载任务在park.MJD2000的状态 (propagate all piggyback missions to the epoch of parking orbit)
                [svCarry , ~] = PropAll_(carry , park.MJD2000 , aux);
                % 停泊轨道状态和所有搭载任务状态在相同时刻求差 (evaluate the state difference between parking orbit and all the piggyback missions at the same epoch)
                temp = repmat(park.sv , size(aux.carry.id , 1) , 1) - svCarry;
                % 计算使用的搭载任务编号 (find the index of piggyback missions)
                tempIndex = find( (sqrt(temp(: , 1).^2 + temp(: , 2).^2 + temp(: , 3).^2) < aux.diffRMax) & (sqrt(temp(: , 4).^2 + temp(: , 5).^2 + temp(: , 6).^2) < aux.diffVMax) );
                
                % 检查: 搭载任务轨道的近地点高度
                if ~isempty(intersect(tempIndex , find(carry.oe(: , 1) .* (1 - carry.oe(: , 2)) - aux.req < aux.altiMin)))
                    if aux.langFlag == 1
                        disp('错误：搭载任务的近地点高度低于500km！')
                    else
                        disp('Warning: the apogee altitude of piggyback mission is lower than 500km!')
                    end
                    aux.errorFlag = 1;
                end
                
                % 检查：停泊轨道是否为某条搭载任务轨道 (check: whether the parking orbit is one of the piggyback missions)
                if isempty(tempIndex)
                    if aux.langFlag == 1
                        disp('错误：停泊轨道非搭载任务轨道！')
                    else
                        disp('Warning: the parking orbit is not one of the piggyback missions!')
                    end
                    aux.errorFlag = 1;
                else
                    % 统计搭载使用次数(count which piggyback mission is used)
                    carryAray(tempIndex) = carryAray(tempIndex) + 1 / size(tempIndex , 1);
                end
        end
    end
    
    % 验证出发状态 (check the state before the first impulse)
    if strcmp(Flag , 'From')
        
        % 子任务数 + 1 (number of subtask + 1)
        NsubTask = NsubTask + 1;
        
        % 提取目标导航卫星编号，编号应和'constellation.txt'中的卫星相对应
        % (extract the index of navigation satellite, and it should be same with the index listed in 'constellation.txt')
        tline = fgets(fid);
        lineTemp = textscan(tline , '%s %s %f');
        tarIndex = lineTemp{3};
        
        % 读取子任务出发时的状态 (read the state before the first impulse)
        tline = fgets(fid);
        tline = fgets(fid);
        xcraft = cell2mat(textscan(tline , '%f'))';
        
        % 检查: 子任务出发时间应在建设区间内，且应晚于停泊轨道时间 (check: whether the epoch of first impulse is located in the construction phase, and later than the epoch of parking orbit)
        if xcraft(1) < aux.MJD2000_ti || xcraft(1) > aux.MJD2000_tf || xcraft(1) < park.MJD2000
            if aux.langFlag == 1
                disp(['错误：导航卫星',num2str(tarIndex),'第一次脉冲时间早于停泊轨道时间或不在星座建设时间内！']);
            else
                disp(['Warning: for satellite ',num2str(tarIndex),', the epoch of first impulse is earlier than the epoch of parking orbit or out of the construction phase!']);
            end
            aux.errorFlag = 1;
        end
        
        % 卫星在停泊轨道出发前状态 (evaluate the state before the first impulse)
        [svTemp , ~] = PropOne_(park , xcraft(1) , aux , 1);
        
        % 计算卫星在第一次脉冲前的位置，速度和质量误差 (evaluate position and velocity difference before the first impulse)
        diffR = norm(xcraft(2:4) - svTemp(1:3));
        diffV = norm(xcraft(5:7) - svTemp(4:6));
        diffM = xcraft(8) - aux.m0;
        
        % 检查: 判断出发状态是否满足精度要求 (check: whether the difference meet the constraints)
        if diffR > aux.diffRMax || diffV > aux.diffVMax || diffM > aux.diffMMax
            if aux.langFlag == 1
                disp(['错误：卫星',num2str(tarIndex),'在第一次脉冲前的状态违反约束!'])
            else
                disp(['Warning：for satellite ',num2str(tarIndex),', the state before the first impulse violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
    end
    
    % 验证脉冲弧段 (check the state of impulse)
    if strcmp(Flag , 'Impulse')
        
        % 读取脉冲状态 (read the state of impulse)
        tline = fgets(fid);
        Impulse = textscan(tline , '%f');
        Impulse = cell2mat(Impulse)';
        
        % 验证轨道连续性 (verify the continuity of state)
        aux = verifyConsist([Impulse(1 : 7) , xcraft(8)] , xcraft , aux);
        
        % 验证变轨质量 (evaluate the mass after impulse)
        dv = Impulse(8 : 10);
        mcraft0 = xcraft(8) ;
        mcraft1 = Impulse(11);
        mf = mcraft0 * exp(-norm(dv) / aux.c);
        
        % 检查：脉冲后卫星质量 (check: whether the mass after impulse meet the constraint)
        if abs(mf - mcraft1) > aux.diffMMax
            if aux.langFlag == 1
                disp(['错误：卫星',num2str(tarIndex),'变轨质量计算错误!'])
            else
                disp(['Warning：for satellite ',num2str(tarIndex),', the mass after impulse violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
        
        % 更新飞行器状态 (update current state)
        xcraft(1:4) = Impulse(1:4);
        xcraft(5:7) = Impulse(5:7) + Impulse(8:10);
        xcraft(8) = mcraft1;
    end
    
    % 验证滑行弧段 (check the coast segment)
    if strcmp(Flag , 'Coast')
        
        % 读取滑行段的状态 (read the state of coast segment)
        tline = fgets(fid);
        state0 = textscan(tline , '%f'); % 滑行开始时卫星的状态 (the state of coast segment begins)
        state0 = cell2mat(state0)';
        tline = fgets(fid);
        state1 = textscan(tline , '%f'); % 滑行结束时卫星的状态 (the state of coast segment ends)
        state1 = cell2mat(state1)';
        
        % 验证轨道连续性 (verify the continuity of state)
        aux = verifyConsist(state0 , xcraft , aux);
        
        % 验证滑行段积分精度 (evaluate the positon and velocity error with numerical integration)
        [diffR , diffV] = verifyCoast(state0 , state1 , aux);
        
        % 检查：滑行段轨道积分误差 (check: whether the integration error meet the constraint)
        if diffR > aux.diffRMax || diffV > aux.diffVMax
            if aux.langFlag == 1
                disp(['错误：卫星',num2str(tarIndex),'滑行段轨道积分误差违反约束!'])
            else
                disp(['Warning：for satellite ',num2str(tarIndex),', the propagation error of coast segment violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
        
        % 更新飞行器状态 (update current state)
        xcraft = state1;
    end
    
    % 验证交会状态 (check the state of injection)
    if strcmp(Flag , 'Injection')
        
        % 读取入轨状态 (read the state of injection)
        tline = fgets(fid);
        Inject = textscan(tline , '%f');
        Inject = cell2mat(Inject)';
        
        % 计算交会时刻导航卫星状态 (evaluate the state of satellite after injection)
        [svTar , ~] = PropOne_(cons , Inject(1) , aux , tarIndex);
        
        % 验证状态连续性 (verify the continuity of state)
        aux = verifyConsist(Inject , xcraft , aux);
        
        % 验证入轨时刻的位置与速度偏差 (evaluate the position and velocity difference after injection)
        diffR = norm(Inject(2 : 4) - svTar(1 : 3));
        diffV = norm(Inject(5 : 7) - svTar(4 : 6));
        
        % 检查：入轨位置和速度误差 (check: whether the position and velocity difference meet the constraint)
        if diffR > aux.diffRMax || diffV > aux.diffVMax
            if aux.langFlag == 1
                disp(['错误：卫星',num2str(tarIndex),'入轨精度不满足!'])
            else
                disp(['Warning：for satellite ',num2str(tarIndex),', the injection error violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
        
        % 检查：入轨后剩余质量约束 (check: whether the mass after injection greater than dry mass)
        if Inject(8) < aux.mdry
            if aux.langFlag == 1
                disp(['错误：卫星',num2str(tarIndex),'入轨后质量小于干重!'])
            else
                disp(['Warning：for satellite ',num2str(tarIndex),', the residual mass is lower than the dry mass!'])
            end
            aux.errorFlag = 1;
        end
        
        % 显示验证结果 (display verification result of injection)
        if aux.errorFlag == 0
            if aux.langFlag == 1
                disp(['正确：导航卫星', num2str(tarIndex),'成功入轨！'])
            else
                disp(['Correct：satellite ',num2str(tarIndex),' is injected into the orbit!'])
            end
        end
        
        % 成功入轨后visitAray相应位置清零 (if the satellite is injected into the target orbit successfully,
        % the corresponding element in the visitAray is replaced by 0)
        if aux.errorFlag == 0
            visitAray(tarIndex) = 0;
        else
            aux.errorFlag = 0;
        end
        
        % 检查: 两种部署方式的导航卫星数是否满足约束 (check: whether the number of satellites for two kind of deployments meet the constraint)
        switch typeFlag
            case 1 % 发射新火箭 (for new launchings)
                if NsubTask > 16
                    if aux.langFlag == 1
                        disp('错误：发射新火箭方式的导航卫星数大于16！');
                    else
                        disp('Warning: for launching new rocket, the navigation satellites are greater than 16!');
                    end
                    aux.errorFlag = 1;
                end
            case 2 % 利用搭载 (for piggyback missions)
                if NsubTask > 8
                    if aux.langFlag == 1
                        disp('错误：搭载方式的导航卫星数大于8！');
                    else
                        disp('Warning: for piggyback mission, the navigation satellites are greater than 8!');
                    end
                    aux.errorFlag = 1;
                end
        end
        
        % 入轨导航卫星数 + 1 (number of injected satellites + 1)
        Nnavi = Nnavi + 1;
        
    end
    
    % 读取一行 (read a new line)
    tline = fgets(fid);
    
end

% 检查: 如果搭载任务被多次使用 (check: whether piggyback missions are used more than once)
if any(carryAray > 1)
    if aux.langFlag == 1
        disp('----------------------------------------------------------');
        disp(['错误：搭载任务被多次使用！编号为: ' , num2str(find(carryAray > 1)')])
    else
        disp('----------------------------------------------------------');
        disp(['Warning：the piggyback missions are used more than once！indexes: ' , num2str(find(visitAray > 1)')])
    end
    aux.errorFlag = 1;
end

% 检查: 是否还有未入轨卫星 (check: whether some satellites failded to inject into the target orbits)
if any(visitAray)
    if aux.langFlag == 1
        disp('----------------------------------------------------------');
        disp(['错误：部分卫星未能正确入轨！编号为: ' , num2str(find(visitAray ~= 0)')])
    else
        disp('----------------------------------------------------------');
        disp(['Warning：some satellites failed to inject into the target orbits！indexes: ' , num2str(find(visitAray ~= 0)')])
    end
    aux.errorFlag = 1;
end

% 计算Obj2 (evaluate Obj2)
Obj2 = Nlaunch * aux.cLaunch + Ncarry * aux.cCarry + Nnavi * aux.cNavi;
if Obj2 > aux.minCost
    if aux.langFlag == 1
        disp(['错误：星座建设成本' , num2str(Obj2) , '大于约束!']);
    else
        disp(['Warning：the cost of constellation violates the constraint!']);
    end
    aux.errorFlag = 1;
end

% 显示'transfer.txt'的验证结果 (display the verification result of 'transfer.txt')
if aux.errorFlag == 0
    if aux.langFlag == 1
        disp('transfer.txt 通过验证！')
    else
        disp('transfer.txt passed the verification！')
    end
end

% 输出结果 (output verification results)
if aux.errorFlag == 0
    if aux.langFlag == 1
        disp('==========================================================');
        disp('通过所有验证，输出星座设计结果：')
        disp(['Obj1: ' , num2str(Obj1)])
        disp(['Obj2: ' , num2str(Obj2) , ' (货币单位)'])
        disp(['发射火箭数: ' , num2str(Nlaunch)])
        disp(['使用搭载数: ' , num2str(Ncarry)])
        disp(['导航卫星数: ' , num2str(Nnavi)])
    else
        disp('==========================================================');
        disp('The results are listed below：')
        disp(['Obj1: ' , num2str(Obj1)])
        disp(['Obj2: ' , num2str(Obj2) , ' (currency unit)'])
        disp(['Number of new launches: ' , num2str(Nlaunch)])
        disp(['Number of piggybacks: ' , num2str(Ncarry)])
        disp(['Number of satellites: ' , num2str(Nnavi)])
    end
else
    if aux.langFlag == 1
        disp('==========================================================');
        disp('未能通过所有验证，请检查相关错误！')
    else
        disp('==========================================================');
        disp('Failed to pass all validations，please check errors！')
    end
end

end


% =========================================================================
% 构造仿真参数结构体
%
% 输出：
% aux              [structure]          仿真参数结构体
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Generate a structure of simulation parameters
%
% Output：
% aux              [structure]          structure of simulation parameters
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function aux = aux_
aux.mu = 398600; % 地球引力常数(km^3/s^2) (Gravitational constant of the Earth, km^3/s^2)
aux.req = 6378; % 地球半径(km) (Equatorial radius, km)
aux.j2 = 0.0010826; % j2摄动(--) (Earth oblateness gravity coefficient, --)
aux.we = 7.29211585530e-5; % 地球自转速率(rad/sec) (Earth rotation rate, rad/sec)
aux.g0 = 9.80665; % 海平面重力加速度(m/s^2) (Gravity acceleration at sea level, m/s^2)
aux.Isp = 300; % 比冲(s) (Specific impulse, s)
aux.c = 1e-3 * aux.g0 * aux.Isp; % 燃气速度 (Gas velocity, km/s)
aux.m0 = 100; % 卫星初始重量(kg) (Initial mass of navigation satellite, kg)
aux.mdry = 80; % 卫星干重(kg) (Dry mass of navigation satellite, kg)
aux.eleMin = 10 * pi / 180; % 城市最低可视仰角(rad) (Minimum elevation angle, rad)
aux.dtr = pi / 180; % 度转弧度 (deg to rad)
aux.rtd = 180 / pi; % 弧度转度 (rad to deg)
aux.RelTol = 1e-13; % 数值积分绝对误差 (Relative integration tolerance)
aux.AbsTol = 1e-13; % 数值积分相对误差 (Absolute integration tolerance)
aux.diffTMax = 1e-3; % 验证程序时间精度(s) (time tolerance, s)
aux.diffRMax = 1e-3; % 验证程序位置精度(km) (position tolerance, km)
aux.diffVMax = 1e-6; % 验证程序速度精度(km/s) (velocity tolerance, km/s)
aux.diffMMax = 1e-4; % 验证程序质量精度(kg) (mass tolerance, kg)
aux.tStep = 120; % 服务阶段采样步长(s) (Sampling steps in service stage, s)
aux.JDtoMJD2000 = -2451544.5; % 从JD到MJD2000(day) (from JD to MJD2000, day)
aux.MJD2000toJD = 2451544.5; % 从MJD2000到JD(day) (from MJD2000 to JD, day)
aux.MJD2000_ti = jday(2020 , 1 , 1 , 0 , 0 , 0) + aux.JDtoMJD2000; % 建设阶段开始修正儒略日(day) (epoch of construction begins, day)
aux.MJD2000_tf = jday(2020 , 4 , 1 , 0 , 0 , 0) + aux.JDtoMJD2000; % 建设阶段结束修正儒略日(day) (epoch of construction ends, day)
aux.cLaunch = 1.2; % 单次火箭成本(货币单位) (Single launch cost, currency unit)
aux.cCarry = 0.2; % 单次搭载成本(货币单位) (Single piggyback cost, currency unit)
aux.cNavi = 0.05; % 单颗导航卫星成本(货币单位) (Single navigation satellite cost, currency unit)
aux.minCost = 10; % 星座建设成本(货币单位) (Constellation construction cost, currency unit)
aux.maxGDOP = 10; % 城市GDOP约束 (GDOP constraint for each city)
aux.altiMin = 500; % 轨道高度约束 (Minimum altitude constraint, km)
aux.altiPark = 900; % 火箭入轨高度 (Altitude of parking orbit for launching new rockets, km)
aux.errorFlag = 0; % 错误标示 (Error flag during verification)
tt1 = [86400 * 0 : aux.tStep : 86400 * 1]; % 第1天(Sampling points of 1-th day, s)
tt2 = [86400 * 6 : aux.tStep : 86400 * 7]; % 第7天(Sampling points of 7-th day, s)
tt3 = [86400 * 29 : aux.tStep : 86400 * 30]; % 第30天(Sampling points of 30-th day, s)
aux.tAray = [tt1 , tt2 , tt3]; % 服务阶段开始后的（即从MJD2000=7396开始）采样时间点, (Sampling points from MJD2000 = 7396, s)
end

% =========================================================================
% 计算GDOP
%
% 输入:
% satMtx          [nx3]          可见卫星在ecf系下的位置坐标
% obsEcf          [1x3]          城市在ecf系下的位置坐标
%
% 输出:
% gdop            [1x1]           几何精度因子(GDOP)
%
% 测试：
% satMtx = [15524.471175 , -16649.826222 , 13512.2723887;
%          -2304.058534 , -23287.906465 , 11917.038105;
%           16680.243357 , -3069.625561 , 20378.551047;
%          -14799.931395 , -21425.35824 , 6069.947224];
% obsEcf = [-730.000 , -5440.000 , 3230.000];
% gdop_(satMtx , obsEcf)
%
% 参考文献：
% https://en.wikipedia.org/wiki/Dilution_of_precision_(navigation)
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Evaluate GDOP
%
% Input:
% satMtx          [nx3]          satellite positions under ECF frame
% obsEcf          [1x3]          city position under ECF frame
%
% Output:
% gdop            [1x1]          Geometric Dilution Precision (GDOP)
%
% Example：
% satMtx = [15524.471175 , -16649.826222 , 13512.2723887;
%     -2304.058534 , -23287.906465 , 11917.038105;
%     16680.243357 , -3069.625561 , 20378.551047;
%     -14799.931395 , -21425.35824 , 6069.947224];
% obsEcf = [-730.000 , -5440.000 , 3230.000];
% gdop_(satMtx , obsEcf)
%
% Reference：
% https://en.wikipedia.org/wiki/Dilution_of_precision_(navigation)
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function gdop = getGDOP(satMtx , obsEcf)
seeNum = size(satMtx , 1); % 可见卫星数 (number of visible satellites)
rTemp = satMtx - repmat(obsEcf , seeNum , 1); % 地面站指向卫星的向量 (vector from city to satellites)
rNorm = rTemp ./ repmat(sqrt(rTemp(: , 1).^2 + rTemp(: , 2).^2 + rTemp(: , 3).^2) , 1 , 3); % 地面站指向卫星的单位向量 (unit vector from city to satellites)
H = [rNorm , ones(seeNum , 1)]; % 构造H矩阵 (generate H matrix)
gdop = sqrt( trace(pinv(H' * H)) ); % 计算gdop (evaluate gdop)
end

% =========================================================================
% 计算指标1
%
% 输入:
% aux                  [structure]         仿真参数
%    aux.MJD2000_tf    [1x1]               建设阶段结束修正儒略日(day)
%    aux.tAray         [1x2163]            服务阶段开始后的（即从MJD2000=7396开始）采样时间点(s)
%    aux.city          [structure]         城市信息结构体
%    aux.cons          [structure]         星座信息结构体
%    aux.eleMin        [1x1]               城市最低可视仰角(rad)
%    aux.maxGDOP       [1x1]               城市GDOP约束
%
% 输出:
% Obj1                 [1x1]               设计指标1
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Evaluate Obj1
%
% Input:
% aux                  [structure]         simulation parameters
%    aux.MJD2000_tf    [1x1]               epoch of construction ends, day
%    aux.tAray         [1x2163]            sampling points from MJD2000 = 7396, s
%    aux.city          [structure]         structure of city information
%    aux.cons          [structure]         structure of constellation information
%    aux.eleMin        [1x1]               minimum elevation angle, rad
%    aux.maxGDOP       [1x1]               maximum GDOP constraint for each city
%
% Output:
% Obj1                 [1x1]               the first index
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function Obj1 = getObj1_(aux)
Nt = size(aux.tAray , 2); % 时间采样点数 (number of sampling points)
Ncity = size(aux.city.ecf , 1); % 城市数 (number of cities)
NnaviSat = size(aux.cons.oe , 1); % 导航卫星数 (number of navigation satellites)
eleMin = aux.eleMin; % 城市最低可见仰角 (minimum elevation angle)

% 计算星座在ecf系下的位置速度 (evaluate position and velocity under ECF frame)
consEcfMtx = zeros(NnaviSat * 6 , Nt);
count = 1;
for i = 1 : NnaviSat % 卫星循环 (loop of navigation satellites)
    % 计算ecf下位置速度 (evaluate position and velocity under ECF frame)
    [rv_eci , rv_ecf] = oe2eciecf_(aux.cons.oe(i , :) , aux.MJD2000_tf , aux.MJD2000_tf + aux.tAray / 86400 , aux);
    % 保存星座位置速度 (save position and velocity under ECF frame)
    consEcfMtx( (count - 1) * 6 + 1 : count * 6 , : ) = rv_ecf;
    count = count + 1;
end

% 计算所有城市在所有采样点上的GDOP值 (evaluate GDOP matrix for each city and each sampling point)
GDOPMtx = 1e8 * ones(Ncity , Nt); % 初始化GDOP矩阵 (pre allocate GDOP matrix)

% 时间循环 (time loop)
for j = 1 : Nt
    
    % 城市循环 (city loop)
    for k = 1 : Ncity
        
        % 提取第k个城市坐标, ECF系 (extract position of the k-th city, ECF frame)
        cityEcf = aux.city.ecf(k , :);
        
        % 提取星座所有卫星在第j个采样点的位置速度，ECF系 (extract position of constellation at the j-th sampling point, ECF frame)
        satEcf = reshape(consEcfMtx(: , j) , 6 , NnaviSat)';
        
        % 计算从第k个城市到星座所有卫星的位置矢量 (evaluate vectors from the k-th city to all satellites of constellation)
        city2sat = [satEcf(: , 1) - cityEcf(1) , satEcf(: , 2) - cityEcf(2) , satEcf(: , 3) - cityEcf(3)];
        
        % 计算卫星相对城市的夹角 (evalute angles between cityEcf vector and city2sat vectors)
        thetaAray = acos( (cityEcf(1) .* city2sat(: , 1) + cityEcf(2) .* city2sat(: , 2) + cityEcf(3) .* city2sat(: , 3)) ./ ...
            ( sqrt(cityEcf(1)^2 + cityEcf(2)^2 + cityEcf(3)^2) * sqrt(city2sat(: , 1).^2 + city2sat(: , 2).^2 + city2sat(: , 3).^2) ) );
        
        % 提取可见卫星编号 (extract visible satellite indexes)
        seeIndex = find(thetaAray < (0.5 * pi - eleMin));
        
        % 如果可见卫星数 >= 4，计算GDOP (if the number of visiable satellties is greater than 4, evaluate GDOP)
        if size(seeIndex , 1) >= 4
            GDOPMtx(k , j) = getGDOP(satEcf(seeIndex , 1 : 3) , cityEcf);
        end
        
    end
    
end

% 计算目标函数 (evaluate Obj1)
Obj1 = ( max(GDOPMtx') <= aux.maxGDOP ) * aux.city.weight;

end


% =========================================================================
% 通过偏心率和平近点角计算偏近点角和真近点角
%
% 输入:
% e            [1x1]         偏心率
% M            [nx1]         平近点角
%
% 输出:
% Ek1          [nx1]         偏近点角
% theta        [nx1]         真近点角
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Evaluate eccentric anomaly and true anomaly with eccentricity and mean
% anomaly
%
% Input:
% e            [1x1]         eccentricity
% M            [nx1]         mean anomaly
%
% Output:
% Ek1          [nx1]         eccentric anomaly
% theta        [nx1]         true anomaly
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function [E1 , tanom] = m2ef_(ecc , manom)
error = 1e-12;
E = manom;
E1 = E - (E - ecc .* sin(E) - manom) ./ (1 - ecc .* cos(E));
while any(abs(E1 - E) > error)
    E = E1;
    E1 = E - (E - ecc .* sin(E) - manom) ./ (1 - ecc .* cos(E));
end
tanom = 2 * atan2(sqrt(1 + ecc) .* sin(E1 / 2) , sqrt(1 - ecc) .* cos(E1 / 2));
end

% =========================================================================
% J2摄动轨道动力学模型
%
% 输入:
% t             [1x1]            时间
% x             [6x1]            状态
% aux           [structure]      仿真参数
%   aux.mu      [1x1]            地球引力常数(km^3/s^2)
%   aux.req     [1x1]            地球半径(km)
%   aux.j2      [1x1]            j2摄动(--)
%
% 输出:
% dx            [6x1]            向量场
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Dynamic model under j2 perturbation
%
% Input：
% t             [1x1]             time
% x             [6x1]             state
% aux           [structure]       simulation parameters
%   aux.mu      [1x1]             Gravitational constant of the Earth (km^3/s^2)
%   aux.req     [1x1]             Equatorial radius (km)
%   aux.j2      [1x1]             Earth oblateness gravity coefficient (--)
%
% Output:
% dx            [6x1]             vector field
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function dx = modelJ2(t , x , aux)
mu = aux.mu;
j2 = aux.j2;
req = aux.req;
r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
dx = [x(4);
    x(5);
    x(6);
    -mu * x(1) / r^3 * (1 + 3/2 * j2 * (req/r)^2 * (1 - 5 * x(3)^2 / r^2));
    -mu * x(2) / r^3 * (1 + 3/2 * j2 * (req/r)^2 * (1 - 5 * x(3)^2 / r^2));
    -mu * x(3) / r^3 * (1 + 3/2 * j2 * (req/r)^2 * (3 - 5 * x(3)^2 / r^2))];
end

% =========================================================================
% 轨道根数转eci和ecf下的位置速度
%
% 输入:
% xx_oe             [1x6]           t0时刻的轨道根数
%                                   [sma(km) , ecc(--) , incl(rad) , ...
%                                    raan(rad) , argper(rad) , manom(rad)]
% t0                [1x1]           xx_oe的历元, MJD2000 (day)
% tt                [1xn]           需要计算的时间点, MJD2000 (day)
% aux               [structure]     仿真参数
%   aux.mu          [1x1]           地球引力常数(km^3/s^2)
%   aux.req         [1x1]           地球半径(km)
%   aux.j2          [1x1]           j2摄动(--)
%   aux.we          [1x1]           地球自转速率(rad/sec)
%   aux.MJD2000toJD [1x1]           从MJD2000到JD(day)
%
% 输出：
% rv_eci            [6xn]           tt时刻eci位置速度 [x,y,z,vx,vy,vz] (km, km/s)
% rv_ecf            [6xn]           tt时刻ecf位置速度 [x,y,z,vx,vy,vz] (km, km/s)
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01 (Sep/01/2017)
% -------------------------------------------------------------------------
% Orbital elements to position and velocity under both ECI and ECF frame
%
% Input：
% xx_oe             [1x6]           orbital elements at t0
%                                   [sma(km) , ecc(--) , incl(rad) , ...
%                                    raan(rad) , argper(rad) , manom(rad)]
% t0                [1x1]           epoch of xx_oe, MJD2000 (day)
% tt                [1xn]           sampling points, MJD2000 (day)
% aux               [structure]     simulation parameters
%   aux.mu          [1x1]           Gravitational constant of the Earth (km^3/s^2)
%   aux.req         [1x1]           Equatorial radius (km)
%   aux.j2          [1x1]           Earth oblateness gravity coefficient (--)
%   aux.we          [1x1]           Earth rotation rate (rad/sec)
%   aux.MJD2000toJD [1x1]           from MJD2000 to JD (day)
%
% Output：
% rv_eci            [6xn]           states under ECI frame at sampling points tt [x,y,z,vx,vy,vz] (km, km/s)
% rv_ecf            [6xn]           states under ECF frame at sampling points tt [x,y,z,vx,vy,vz] (km, km/s)
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function [rv_eci , rv_ecf] = oe2eciecf_(xx_oe , t0 , tt , aux)

% 读取参数 (extract simulation parameters)
mu = aux.mu;
j2 = aux.j2;
req = aux.req;
we = aux.we;

% 提取t0时刻的轨道根数 (extract orbital elements at t0)
sma = xx_oe(1); % 半长轴 (km) (semi-major axis, km)
ecc = xx_oe(2); % 偏心率 (--) (eccentricity, --)
incl = xx_oe(3); % 轨道倾角 (rad) (inclination, rad)
raan0 = xx_oe(4); % 升交点赤经(rad) (right ascension of ascending node, rad)
argper0 = xx_oe(5); % 近地点幅角 (rad) (argument of perigee, rad)
manom0 = xx_oe(6); % 平近点角 (rad) (mean anomaly, rad)

tt = tt(:); % 需要计算的时间tt变为列向量 (transform tt to column vector)

% 计算j2产生的平根变化率 (evaluate the change rate of elements under j2 perturbation)
t1 = 1 - ecc.^2;
t2 = (req ./ sma).^2;
An = sqrt(mu ./ sma.^3);
draan = -1.5 * An * j2 ./ t1.^2 .* t2 .* cos(incl);
dargper = 1.5 * An * j2 ./ t1.^2 .* t2 .* (2 - 2.5 * sin(incl).^2);
dmanom = 0.75 * An * j2 ./ t1.^(1.5) .* t2 .* (3 * cos(incl).^2 - 1);

% 计算采样点处的平跟 (evaluate the orbital elements at sampling points)
raan = mod(raan0 + draan .* (tt - t0) * 86400 , 2 * pi);
argper = mod(argper0 + dargper .* (tt - t0)  * 86400 , 2 * pi);
manom = mod(manom0 + (An + dmanom) .* (tt - t0)  * 86400 , 2 * pi);

% 计算真近点角 (evaluate true anomaly)
[~ , tanom] = m2ef_(ecc , manom);

% 计算eci位置速度 (evaluate position and velocity under ECI frame)
r = sma .* (1 - ecc .^ 2) ./ (1 + ecc .* cos(tanom));
v = sqrt( 2 * mu ./ r - mu ./ sma );
gamma = atan( ecc .* sin(tanom) ./ (1 + ecc .* cos(tanom)) );
x = r .* ( cos(tanom + argper) .* cos(raan) - sin(tanom + argper) .* cos(incl) .* sin(raan) );
y = r .* ( cos(tanom + argper) .* sin(raan) + sin(tanom + argper) .* cos(incl) .* cos(raan) );
z = r .* ( sin(tanom + argper) .* sin(incl) );
vx = v .* ( -sin(tanom + argper - gamma) .* cos(raan) - cos(tanom + argper - gamma) .* cos(incl) .* sin(raan) );
vy = v .* ( -sin(tanom + argper - gamma) .* sin(raan) + cos(tanom + argper - gamma) .* cos(incl) .* cos(raan) );
vz = v .* ( cos(tanom + argper - gamma) .* sin(incl));

% 构造eci下的位置速度 (save position and velocity under ECI frame)
r_eci = [x , y , z]';
v_eci = [vx , vy , vz]';
rv_eci = [r_eci ; v_eci];

% eci转ecf (transform eci to ecf)
if nargout > 1
    % 计算采样点tt的格林威治恒星时(rad) (evaluate Greenwich sidereal time, rad)
    vgast = gstime(tt + aux.MJD2000toJD);
    
    % 计算坐标旋转矩阵 (evaluate coordinate transformation matrix)
    tm = zeros(3 , 3 , length(tt));
    tmdot = zeros(3 , 3 , length(tt));
    tm(1 , 1 , :) = cos(vgast); tm(1 , 2 , :) = sin(vgast); tm(1 , 3 , :) = 0;
    tm(2 , 1 , :) = -sin(vgast); tm(2 , 2 , :) = cos(vgast); tm(2 , 3 , :) = 0;
    tm(3 , 1 , :) = 0; tm(3 , 2 , :) = 0; tm(3 , 3 , :) = 1;
    tmdot(1 , 1 , :) = -we * sin(vgast); tmdot(1 , 2 , :) = we * cos(vgast); tmdot(1 , 3 , :) = 0;
    tmdot(2 , 1 , :) = -we * cos(vgast); tmdot(2 , 2 , :) = -we * sin(vgast); tmdot(2 , 3 , :) = 0;
    tmdot(3 , 1 , :) = 0; tmdot(3 , 2 , :) = 0; tmdot(3 , 3 , :) = 0;
    
    % 构造ecf下的位置速度 (evaluate position and velocity under ECF frame)
    r_ecf = x_multi_y_(tm , r_eci);
    v_ecf = x_multi_y_(tm , v_eci) + x_multi_y_(tmdot , r_eci);
    rv_ecf = [r_ecf ; v_ecf];
end

end

% =========================================================================
% J2平跟递推所有卫星在epochAray的位置速度和轨道根数
%
% 输入:
% sat                 [structure]         卫星参数结构体
%    sat.id           [nx1]               卫星编号
%    sat.MJD2000      [nx1]               卫星历元
%    sat.oe           [nx6]               卫星轨道根数
% epochAray           [1xm]               需要计算的时刻
% aux                 [structure]         仿真参数
%    aux.mu           [1x1]               地球引力常数(km^3/s^2)
%    aux.req          [1x1]               地球半径(km)
%    aux.j2           [1x1]               j2摄动(--)
%
% 输出:
% svCubic             [nx6xm]             所有卫星在epochAray的位置速度
% oeCubic             [nx6xm]             所有卫星在epochAray的轨道根数
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Evaluate position, velocity and orbital elements of all satellites at
% given time epochAray
%
% Input:
% sat                 [structure]         parameters of satellite
%    sat.id           [nx1]               index of satellite
%    sat.MJD2000      [nx1]               epoch of satellite
%    sat.oe           [nx6]               orbital elements of satellite
% epochAray           [1xm]               given time
% aux                 [structure]         simulation parameters
%    aux.mu           [1x1]               Gravitational constant of the Earth (km^3/s^2)
%    aux.req          [1x1]               Equatorial radius (km)
%    aux.j2           [1x1]               Earth oblateness gravity coefficient (--)
%
% Output:
% svCubic             [nx6xm]             position and velocity of all satellites at epochAray
% oeCubic             [nx6xm]             orbital elements of all satellites at epochAray
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function [svCubic , oeCubic] = PropAll_(sat , epochAray , aux)

% 提取仿真参数 (extract simulation parameters)
mu = aux.mu;
j2 = aux.j2;
req = aux.req;

% 提取需要计算的时间数和卫星总数 (extract number of sampling points and number of satellites)
epochNum = size(epochAray , 2);
satNum = size(sat.oe , 1);

% 提取所有卫星的参数 (extract epoch and orbital elements of all satellites)
epoch = sat.MJD2000;
sma = sat.oe(: , 1);
ecc = sat.oe(: , 2);
incl = sat.oe(: , 3);
raan0 = sat.oe(: , 4);
argper0 = sat.oe(: , 5);
manom0 = sat.oe(: , 6);

% 计算j2产生的平根变化率 (evaluate the change rate of elements under j2 perturbation)
t1 = 1 - ecc.^2;
t2 = (req ./ sma).^2;
An = sqrt(mu ./ sma.^3);
draan = -1.5 * An * j2 ./ t1.^2 .* t2 .* cos(incl);
dargper = 1.5 * An * j2 ./ t1.^2 .* t2 .* (2 - 2.5 * sin(incl).^2);
dmanom = 0.75 * An * j2 ./ t1.^(1.5) .* t2 .* (3 * cos(incl).^2 - 1);

svCubic = zeros(satNum , 6 , epochNum); % 位置速度预分配 (pre allocate space for position and velocity)
oeCubic = zeros(satNum , 6 , epochNum); % 轨道根数预分配 (pre allocate space for orbital elements)
for i = 1 : epochNum
    
    % 提取需要计算的时刻 (extract sampling points)
    epochNow = epochAray(i);
    
    % 计算采样点处的平跟 (evaluate orbital elements at sampling points)
    raan = mod(raan0 + draan .* (epochNow - epoch) * 86400 , 2 * pi);
    argper = mod(argper0 + dargper .* (epochNow - epoch) * 86400 , 2 * pi);
    manom = mod(manom0 + (An + dmanom) .* (epochNow - epoch) * 86400 , 2 * pi);
    
    % 计算真近点角 (evaluate true anomaly)
    [~ , tanom] = m2ef_(ecc , manom);
    
    % 计算位置速度 (evaluate position and velocity)
    r = sma .* (1 - ecc.^2) ./ (1 + ecc .* cos(tanom));
    v = sqrt(2 * mu ./ r - mu ./ sma );
    gamma = atan(ecc .* sin(tanom) ./ (1 + ecc .* cos(tanom)));
    x = r .* (cos(tanom + argper) .* cos(raan) - sin(tanom + argper) .* cos(incl) .* sin(raan));
    y = r .* (cos(tanom + argper) .* sin(raan) + sin(tanom + argper) .* cos(incl) .* cos(raan));
    z = r .* (sin(tanom + argper) .* sin(incl) );
    vx = v .* (-sin(tanom + argper - gamma) .* cos(raan) - cos(tanom + argper - gamma) .* cos(incl) .* sin(raan));
    vy = v .* (-sin(tanom + argper - gamma) .* sin(raan) + cos(tanom + argper - gamma) .* cos(incl) .* cos(raan));
    vz = v .* (cos(tanom + argper - gamma) .* sin(incl));
    
    % 保存位置速度，轨道根数 (save position velocity and orbital elements)
    svCubic(: , : , i) = [x , y , z , vx , vy , vz];
    oeCubic(: , : , i) = [sma , ecc , incl , raan , argper , manom];
end

end

% =========================================================================
% J2平跟递推sat结构体中编号为index的卫星，得到目标历元时刻(tarEpoch)的位置
% 速度和轨道根数
%
% 输入:
% sat                  [structure]        卫星参数
%    sat.id            [nx1]              卫星编号
%    sat.MJD2000       [nx1]              卫星初始历元时刻
%    sat.oe            [nx6]              卫星在初始历元的轨道根数
%    sat.sv            [nx6]              卫星在初始历元的位置速度
% tarEpoch             [1x1]              需要计算的历元时刻
% aux                  [structure]        仿真参数
%    aux.mu            [1x1]              地球引力常数(km^3/s^2)
%    aux.req           [1x1]              地球半径(km)
%    aux.j2            [1x1]              j2摄动(--)
% index                [1x1]              需要计算的卫星编号
%
% 输出:
% sv                   [1x6]              卫星在tarEpoch的位置速度
% oe                   [1x6]              卫星在tarEpoch的轨道根数
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Evaluate position, velocity and orbital elements of the 'index-th' satellite
% in the 'sat' structure at target epoch (tarEpoch)
%
% Input:
% sat                  [structure]        parameters of satellites
%    sat.id            [nx1]              index of satellites
%    sat.MJD2000       [nx1]              epoch of satellites
%    sat.oe            [nx6]              orbital elements of satellites
%    sat.sv            [nx6]              position and velocity of satellites
% tarEpoch             [1x1]              target epoch
% aux                  [structure]        simulation parameters
%    aux.mu            [1x1]              Gravitational constant of the Earth (km^3/s^2)
%    aux.req           [1x1]              Equatorial radius (km)
%    aux.j2            [1x1]              Earth oblateness gravity coefficient (--)
% index                [1x1]              satellite index
%
% Output:
% sv                   [nx6xm]            position and velocity of all satellites at epochAray
% oe                   [nx6xm]            orbital elements of all satellites at epochAray
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function [sv , oe] = PropOne_(sat , tarEpoch , aux , index)

% 提取仿真参数 (extract simulation parameters)
mu = aux.mu;
j2 = aux.j2;
req = aux.req;

% 提取sat中编号为index的卫星在初始历元时刻的轨道根数
epoch = sat.MJD2000(index);
sma = sat.oe(index , 1);
ecc = sat.oe(index , 2);
incl = sat.oe(index , 3);
raan0 = sat.oe(index , 4);
argper0 = sat.oe(index , 5);
manom0 = sat.oe(index , 6);

% 提取需要计算的历元时刻 (extract the target epoch)
tarEpoch = tarEpoch(:);
epochNum = size(tarEpoch , 1);

% 计算j2产生的平根变化率 (evaluate the change rate of elements under j2 perturbation)
t1 = 1 - ecc.^2;
t2 = (req ./ sma).^2;
An = sqrt(mu ./ sma.^3);
draan = -1.5 * An * j2 ./ t1.^2 .* t2 .* cos(incl);
dargper = 1.5 * An * j2 ./ t1.^2 .* t2 .* (2 - 2.5 * sin(incl).^2);
dmanom = 0.75 * An * j2 ./ t1.^(1.5) .* t2 .* (3 * cos(incl).^2 - 1);

% 计算tarEpoch处的平跟 (evaluate orbital elements at sampling points)
raan = mod(raan0 + draan .* (tarEpoch - epoch) * 86400 , 2 * pi);
argper = mod(argper0 + dargper .* (tarEpoch - epoch) * 86400 , 2 * pi);
manom = mod(manom0 + (dmanom + An) .* (tarEpoch - epoch) * 86400 , 2 * pi);

% 计算真近点角 (evaluate true anomaly)
[~ , theta] = m2ef_(ecc , manom);

% 计算位置速度 (evaluate position and velocity)
r = sma .* (1 - ecc.^2) ./ (1 + ecc .* cos(theta));
v = sqrt(2 * mu ./ r - mu ./ sma );
gamma = atan(ecc .* sin(theta) ./ (1 + ecc .* cos(theta)));
x = r .* (cos(theta + argper) .* cos(raan) - sin(theta + argper) .* cos(incl) .* sin(raan));
y = r .* (cos(theta + argper) .* sin(raan) + sin(theta + argper) .* cos(incl) .* cos(raan));
z = r .* (sin(theta + argper) .* sin(incl) );
vx = v .* (-sin(theta + argper - gamma) .* cos(raan) - cos(theta + argper - gamma) .* cos(incl) .* sin(raan));
vy = v .* (-sin(theta + argper - gamma) .* sin(raan) + cos(theta + argper - gamma) .* cos(incl) .* cos(raan));
vz = v .* (cos(theta + argper - gamma) .* sin(incl));

% 保存位置速度，轨道根数 (save position velocity and orbital elements)
sv = [x , y , z , vx , vy , vz];
oe = [sma * ones(epochNum , 1) , ecc * ones(epochNum , 1), incl * ones(epochNum , 1) , raan , argper , manom];

end

% -------------------------------------------------------------------------
% 检查滑行段的积分精度
%
% 输入:
% state0                    [1x8]               滑行开始状态
% state1                    [1x8]               滑行结束状态
% aux                       [structure]         仿真参数
%    aux.RelTol             [1x1]               数值积分绝对误差
%    aux.RelTol             [1x1]               数值积分相对误差
%
% 输出:
% diffR                     [1x1]               位置误差
% diffV                     [1x1]               速度误差
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% verify the positon and velocity error of coast segment
%
% Input:
% state0                    [1x8]               state before coast
% state1                    [1x8]               state after coast
% aux                       [structure]         simulation parameters
%    aux.RelTol             [1x1]               Relative integration tolerance
%    aux.RelTol             [1x1]               Absolute integration tolerance
%
% Output:
% diffR                     [1x1]               position error
% diffV                     [1x1]               velocity error
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function  [diffR , diffV] = verifyCoast(state0 , state1 , aux)
tSpan = [0 , state1(1) - state0(1)] * 86400;
options = odeset('RelTol' , aux.RelTol , 'AbsTol' , aux.AbsTol);
[~ , rvTemp] = ode113(@modelJ2 , tSpan , state0(2:7) , options , aux);
diffR = norm(state1(2 : 4) - rvTemp(end , 1 : 3));
diffV = norm(state1(5 : 7) - rvTemp(end , 4 : 6));
end

% -------------------------------------------------------------------------
% 验证轨道是否连续
%
% 输入:
% coast0                    [1x8]                 第一段状态
% xcraft                    [1x8]                 第二段状态
% aux                       [structure]           仿真参数
%   aux.diffTMax            [1x1]                 验证程序时间精度 (s)
%   aux.diffRMax            [1x1]                 验证程序位置精度 (km)
%   aux.diffVmax            [1x1]                 验证程序速度精度 (km/s)
%   aux.diffMMax            [1x1]                 验证程序质量精度 (kg)
%   aux.langFlag            [1x1]                 语言标识
%   aux.errorFlag           [1x1]                 错误标识
%
% 输出:
% aux                       [structure]           仿真参数
%
% Copyright 张晨(Chen Zhang)
% 清华大学宇航中心(Tsinghua Space Center)
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% check the continuity of two segments
%
% Input:
% coast0                    [1x8]                 state of first segment
% xcraft                    [1x8]                 state of second segment
% aux                       [structure]           simulation parameters
%   aux.diffTMax            [1x1]                 time tolerance (s)
%   aux.diffRMax            [1x1]                 position tolerance (km)
%   aux.diffVmax            [1x1]                 velocity tolerance (km/s)
%   aux.diffMMax            [1x1]                 mass tolerance (kg)
%   aux.langFlag            [1x1]                 language flag
%   aux.errorFlag           [1x1]                 error flag
%
% Output:
% aux                       [structure]           simulation parameters
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function aux = verifyConsist(coast0 , xcraft , aux)
diffTemp = coast0 - xcraft;
dT = abs(diffTemp(1)) * 86400;
dR = norm(diffTemp(2:4));
dV = norm(diffTemp(5:7));
dM = norm(diffTemp(8));
if abs(dT) > aux.diffTMax || dR > aux.diffRMax  || dV > aux.diffVMax  || dM > aux.diffMMax
    if aux.langFlag == 1
        disp('错误：导航卫星的状态不连续！')
    else
        disp(['Warning: the states are not continuous!']);
    end
    aux.errorFlag = 1;
end
end

% =========================================================================
% 用矩阵运算进行坐标转换
%
% 输入:
% x           [3x3xn]       坐标转换矩阵
% y           [3xn]         位置或速度矩阵
%
% 输出:
% x_multi_y   [3xn]         坐标转换结果
%
% 作者: 张晨, 清华大学宇航中心
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Perform coordinate transformation with matrix operation
%
% Input:
% x           [3x3xn]        coordinate transformation matrix
% y           [3xn]          position or velocity matrix
%
% Output:
% x_multi_y   [3xn]          result matrix after coordinate transformation
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function x_multi_y = x_multi_y_(x , y)
x_multi_y = [sum(permute(x(1,:,:),[3 2 1]) .* y' , 2),...
    sum(permute(x(2,:,:),[3 2 1]) .* y' , 2),...
    sum(permute(x(3,:,:),[3 2 1]) .* y' , 2)]';
end
