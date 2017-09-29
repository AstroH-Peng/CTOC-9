% =========================================================================
% ��9��ȫ���ռ�����ƾ���������Ŀ��֤����
% ��'constellation.txt'��'transfer.txt'���м��
%
% �ο�:
% ������, ��8��ȫ���ռ�����ƾ���������Ŀ��֤����
%
% ����: �ų�, �廪��ѧ�����
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

% ��ʼ������ (initialize parameters)
aux = aux_;

% ѡ����ʾ����(select language)
aux.langFlag = 1; % ����
% aux.langFlag = 2; % English

if aux.langFlag == 1
    disp('==============================================================');
    disp('(1) ���ڼ�� constellation.txt')
else
    disp('==============================================================');
    disp('(1) checking constellation.txt')
end

% ��ȡ city.txt (Read city.txt)
[id , lon , lat , weight] = textread('city.txt' , '%s%f%f%f' , 'headerlines' , 1);
city.id = id;
city.lonlat = [lon , lat];
city.weight = weight;
city.ecf = [aux.req * cosd(lat) .* cosd(lon), ...
    aux.req * cosd(lat) .* sind(lon), ...
    aux.req * sind(lat)];
aux.city = city;

% ��ȡ carry.txt (Read carry.txt)
[id , MJD2000 , smaAray , eccAray , inclAray , raanAray , argperAray , manomAray] = ...
    textread('carry.txt' , '%s%f%f%f%f%f%f%f' , 'headerlines' , 1);
tempMtx = [smaAray , eccAray , inclAray , raanAray , argperAray , manomAray];
for i = 1 : size(id , 1)
    carry.id(i , :) = str2num(cell2mat(id(i)));
end
carry.MJD2000 = MJD2000;
carry.oe = tempMtx;
aux.carry = carry;

% ��ȡ constellation.txt (Read constellation.txt)
[id , MJD2000 , smaAray , eccAray , inclAray , raanAray , argperAray , manomAray] = ...
    textread('constellation.txt' , '%s%f%f%f%f%f%f%f' , 'headerlines' , 0);
tempMtx = [smaAray , eccAray , inclAray , raanAray , argperAray , manomAray];
for i = 1 : size(id , 1)
    cons.id(i , :) = str2num(cell2mat(id(i)));
end
cons.MJD2000 = MJD2000;
cons.oe = tempMtx;
aux.cons = cons;

% ��飺������� (check the indexes of the constellation)
if (cons.id(1) ~= 1) || any(diff(cons.id) - 1)
    if aux.langFlag == 1
        disp('����������һ�����Ǳ�Ų�Ϊ1���߱�Ų�������')
    else
        disp('Warning: the first index is not 1 or the constellation indexes are not continuous!')
    end
    aux.errorFlag = 1;
end

% ��飺������Ԫ (check the epoch of the constellation)
if any(MJD2000 - aux.MJD2000_tf)
    if aux.langFlag == 1
        disp('����������ʼ����ʱ�����')
    else
        disp('Warning: the epoch that constellation begins service is wrong!')
    end
    aux.errorFlag = 1;
end

% ��飺��������߶�Լ�� (check the perigee altitude of the constellation)
if any(smaAray .* (1 - eccAray) - aux.req < aux.altiMin)
    if aux.langFlag == 1
        disp('�����������ص�߶ȵ���500km��')
    else
        disp('Warning: the perigee altitude of the constellation is lower than 500km!')
    end
    aux.errorFlag = 1;
end

% ��飺������Բ������Լ�� (check the inclination of the constellation)
if any( min( abs( [inclAray(eccAray ~= 0) - acos(sqrt(1 / 5)) , inclAray(eccAray ~= 0) - pi + acos(sqrt(1 / 5))] )' ) > 1e-6)
    if aux.langFlag == 1
        disp('������������Բ���û�в����ٽ���ǣ�')
    else
        disp('Waring: the critical inclination is not used in elliptic orbits!')
    end
    aux.errorFlag = 1;
end

% ����Obj1 (evaluate Obj1)
Obj1 = getObj1_(aux);

% ��ʾ����� (display the verification result of constellation.txt)
if aux.errorFlag == 0
    if aux.langFlag == 1
        disp('constellation.txt ͨ����֤��');
    else
        disp('constellation.txt passed the verification!');
    end
else
    if aux.langFlag == 1
        disp('constellation.txt δ��ͨ����֤��������ش���')
    else
        disp('constellation.txt failed to pass the verification! please check errors!')
    end
end

if aux.langFlag == 1
    disp('==============================================================');
    disp('(2) ���ڼ�� transfer.txt')
else
    disp('==============================================================');
    disp('(2) checking transfer.txt')
end

visitAray = ones(size(cons.id)); % ͳ��������ǣ�ÿ���һ��������Ӧλ����Ϊ0 (count the number of injected satellites)
carryAray = zeros(size(carry.id)); % ͳ��ʹ�õĴ������� (count the number of piggyback missions)
NTask = 0; % ͳ�������� (count the number of tasks)
NsubTask = 0; % ͳ���������� (count the number of subtasks)
Nnavi = 0; % ͳ�Ƶ��������� (count the number of navigation satellites)
Nlaunch = 0; % ͳ�ƻ���� (count the number of launching rockets)
Ncarry = 0; % ͳ�ƴ����� (count the number of piggyback missions)

fid = fopen(['transfer.txt'] , 'r'); % ��ȡ 'transfer.txt' (read 'transfer.txt')
tline = fgets(fid); % ��ȡ��һ�� (read the first line)
while tline ~= -1 % ��ȡ'transfer.txt'ֱ�����һ�� (read 'transfer.txt' until the end line)
    
    % ��ȡtline (read tline)
    lineTemp = textscan(tline ,'%s %f');
    Flag = lineTemp{1};
    
    % �����⵽'Task' (if 'task' is detected)
    if strcmp(Flag , 'Task')
        
        % �������� + 1 (number of task + 1)
        NTask = NTask + 1;
        
        % ������������ (clear subtask)
        NsubTask = 0;
        
        % ��ȡ�������� (read tline)
        tline = fgets(fid);
        lineTemp = textscan(tline , '%s %s');
        
        % ʶ���������� (identify the type of task)
        if strcmp(lineTemp{2} , 'Launch')
            Nlaunch = Nlaunch + 1; % �������� + 1 (number of new launchings + 1)
            typeFlag = 1;
        else
            Ncarry = Ncarry + 1; % ������ + 1 (number of piggyback missions + 1)
            typeFlag = 2;
        end
        
        % ��ȡͣ�����״̬ (read state of parking orbit)
        tline = fgets(fid);
        tline = fgets(fid);
        lineTemp = cell2mat(textscan(tline , '%f'))';
        
        % ����ͣ������ṹ�� (construct a structure of parking orbit)
        park.MJD2000 = lineTemp(1); % ͣ�������Ԫ (epoch of parking orbit)
        park.oe = lineTemp(2:7); % ͣ������������ (orbital elements of parking orbit)
        [eciTemp , ~] = oe2eciecf_(park.oe , park.MJD2000 , park.MJD2000 , aux); % �������תλ���ٶ� (from orbital elements to position and velocity)
        park.sv = eciTemp'; % ͣ�����λ���ٶ� (save position and velocity of parking orbit)
        
        % ��飺ͣ�������ʼʱ���Ƿ��ڽ���ʱ��Լ���� (check: whether the epoch of parking orbit is located in the construction phase)
        if park.MJD2000 < aux.MJD2000_ti || park.MJD2000  > aux.MJD2000_tf
            if aux.langFlag == 1
                disp('����ͣ�������ʼʱ�䲻����������ʱ��Լ���ڣ�');
            else
                disp('Warning: the epoch of parking orbit is out of the construction phase!');
            end
            aux.errorFlag = 1;
        end
        
        % ��飺ͣ�����ÿ�����������Ƿ����100kg (check: whether the mass of each satellite at parking orbit is equal to 100kg)
        if lineTemp(8) - aux.m0 ~= 0
            if aux.langFlag == 1
                disp('����ͣ�����������������������100kg��');
            else
                disp('Warning��the initiall mass of each satellite at parking orbit is not equal to 100kg!');
            end
            aux.errorFlag = 1;
        end
        
        % ���ͣ�����״̬ (check the state of parking orbit)
        switch typeFlag
            case 1 % �����»�� (for new launchings)
                % ��飺ͣ������߶��Ƿ�Ϊ900km (check: whether the altitude of parking orbit is equal to 900km)
                if abs(park.oe(1) - (aux.altiPark + aux.req)) > aux.diffRMax
                    if aux.langFlag == 1
                        disp('����ͣ������߶�Υ��Լ����')
                    else
                        disp('Warning: the altitude of parking orbit is not equal to 900km!')
                    end
                    aux.errorFlag = 1;
                end
                % ��飺ͣ������Ƿ�ΪԲ��� (check: whether the parking orbit is circular)
                if park.oe(2) > 1e-6
                    if aux.langFlag == 1
                        disp('����ͣ�������Բ�����')
                    else
                        disp('Warning: the parking orbit is not circular!')
                    end
                    aux.errorFlag = 1;
                end
                
            case 2 % ���� (for piggyback missions)
                % �������д���������park.MJD2000��״̬ (propagate all piggyback missions to the epoch of parking orbit)
                [svCarry , ~] = PropAll_(carry , park.MJD2000 , aux);
                % ͣ�����״̬�����д�������״̬����ͬʱ����� (evaluate the state difference between parking orbit and all the piggyback missions at the same epoch)
                temp = repmat(park.sv , size(aux.carry.id , 1) , 1) - svCarry;
                % ����ʹ�õĴ��������� (find the index of piggyback missions)
                tempIndex = find( (sqrt(temp(: , 1).^2 + temp(: , 2).^2 + temp(: , 3).^2) < aux.diffRMax) & (sqrt(temp(: , 4).^2 + temp(: , 5).^2 + temp(: , 6).^2) < aux.diffVMax) );
                
                % ���: �����������Ľ��ص�߶�
                if ~isempty(intersect(tempIndex , find(carry.oe(: , 1) .* (1 - carry.oe(: , 2)) - aux.req < aux.altiMin)))
                    if aux.langFlag == 1
                        disp('���󣺴�������Ľ��ص�߶ȵ���500km��')
                    else
                        disp('Warning: the apogee altitude of piggyback mission is lower than 500km!')
                    end
                    aux.errorFlag = 1;
                end
                
                % ��飺ͣ������Ƿ�Ϊĳ������������ (check: whether the parking orbit is one of the piggyback missions)
                if isempty(tempIndex)
                    if aux.langFlag == 1
                        disp('����ͣ������Ǵ�����������')
                    else
                        disp('Warning: the parking orbit is not one of the piggyback missions!')
                    end
                    aux.errorFlag = 1;
                else
                    % ͳ�ƴ���ʹ�ô���(count which piggyback mission is used)
                    carryAray(tempIndex) = carryAray(tempIndex) + 1 / size(tempIndex , 1);
                end
        end
    end
    
    % ��֤����״̬ (check the state before the first impulse)
    if strcmp(Flag , 'From')
        
        % �������� + 1 (number of subtask + 1)
        NsubTask = NsubTask + 1;
        
        % ��ȡĿ�굼�����Ǳ�ţ����Ӧ��'constellation.txt'�е��������Ӧ
        % (extract the index of navigation satellite, and it should be same with the index listed in 'constellation.txt')
        tline = fgets(fid);
        lineTemp = textscan(tline , '%s %s %f');
        tarIndex = lineTemp{3};
        
        % ��ȡ���������ʱ��״̬ (read the state before the first impulse)
        tline = fgets(fid);
        tline = fgets(fid);
        xcraft = cell2mat(textscan(tline , '%f'))';
        
        % ���: ���������ʱ��Ӧ�ڽ��������ڣ���Ӧ����ͣ�����ʱ�� (check: whether the epoch of first impulse is located in the construction phase, and later than the epoch of parking orbit)
        if xcraft(1) < aux.MJD2000_ti || xcraft(1) > aux.MJD2000_tf || xcraft(1) < park.MJD2000
            if aux.langFlag == 1
                disp(['���󣺵�������',num2str(tarIndex),'��һ������ʱ������ͣ�����ʱ�������������ʱ���ڣ�']);
            else
                disp(['Warning: for satellite ',num2str(tarIndex),', the epoch of first impulse is earlier than the epoch of parking orbit or out of the construction phase!']);
            end
            aux.errorFlag = 1;
        end
        
        % ������ͣ���������ǰ״̬ (evaluate the state before the first impulse)
        [svTemp , ~] = PropOne_(park , xcraft(1) , aux , 1);
        
        % ���������ڵ�һ������ǰ��λ�ã��ٶȺ�������� (evaluate position and velocity difference before the first impulse)
        diffR = norm(xcraft(2:4) - svTemp(1:3));
        diffV = norm(xcraft(5:7) - svTemp(4:6));
        diffM = xcraft(8) - aux.m0;
        
        % ���: �жϳ���״̬�Ƿ����㾫��Ҫ�� (check: whether the difference meet the constraints)
        if diffR > aux.diffRMax || diffV > aux.diffVMax || diffM > aux.diffMMax
            if aux.langFlag == 1
                disp(['��������',num2str(tarIndex),'�ڵ�һ������ǰ��״̬Υ��Լ��!'])
            else
                disp(['Warning��for satellite ',num2str(tarIndex),', the state before the first impulse violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
    end
    
    % ��֤���廡�� (check the state of impulse)
    if strcmp(Flag , 'Impulse')
        
        % ��ȡ����״̬ (read the state of impulse)
        tline = fgets(fid);
        Impulse = textscan(tline , '%f');
        Impulse = cell2mat(Impulse)';
        
        % ��֤��������� (verify the continuity of state)
        aux = verifyConsist([Impulse(1 : 7) , xcraft(8)] , xcraft , aux);
        
        % ��֤������� (evaluate the mass after impulse)
        dv = Impulse(8 : 10);
        mcraft0 = xcraft(8) ;
        mcraft1 = Impulse(11);
        mf = mcraft0 * exp(-norm(dv) / aux.c);
        
        % ��飺������������� (check: whether the mass after impulse meet the constraint)
        if abs(mf - mcraft1) > aux.diffMMax
            if aux.langFlag == 1
                disp(['��������',num2str(tarIndex),'��������������!'])
            else
                disp(['Warning��for satellite ',num2str(tarIndex),', the mass after impulse violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
        
        % ���·�����״̬ (update current state)
        xcraft(1:4) = Impulse(1:4);
        xcraft(5:7) = Impulse(5:7) + Impulse(8:10);
        xcraft(8) = mcraft1;
    end
    
    % ��֤���л��� (check the coast segment)
    if strcmp(Flag , 'Coast')
        
        % ��ȡ���жε�״̬ (read the state of coast segment)
        tline = fgets(fid);
        state0 = textscan(tline , '%f'); % ���п�ʼʱ���ǵ�״̬ (the state of coast segment begins)
        state0 = cell2mat(state0)';
        tline = fgets(fid);
        state1 = textscan(tline , '%f'); % ���н���ʱ���ǵ�״̬ (the state of coast segment ends)
        state1 = cell2mat(state1)';
        
        % ��֤��������� (verify the continuity of state)
        aux = verifyConsist(state0 , xcraft , aux);
        
        % ��֤���жλ��־��� (evaluate the positon and velocity error with numerical integration)
        [diffR , diffV] = verifyCoast(state0 , state1 , aux);
        
        % ��飺���жι��������� (check: whether the integration error meet the constraint)
        if diffR > aux.diffRMax || diffV > aux.diffVMax
            if aux.langFlag == 1
                disp(['��������',num2str(tarIndex),'���жι���������Υ��Լ��!'])
            else
                disp(['Warning��for satellite ',num2str(tarIndex),', the propagation error of coast segment violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
        
        % ���·�����״̬ (update current state)
        xcraft = state1;
    end
    
    % ��֤����״̬ (check the state of injection)
    if strcmp(Flag , 'Injection')
        
        % ��ȡ���״̬ (read the state of injection)
        tline = fgets(fid);
        Inject = textscan(tline , '%f');
        Inject = cell2mat(Inject)';
        
        % ���㽻��ʱ�̵�������״̬ (evaluate the state of satellite after injection)
        [svTar , ~] = PropOne_(cons , Inject(1) , aux , tarIndex);
        
        % ��֤״̬������ (verify the continuity of state)
        aux = verifyConsist(Inject , xcraft , aux);
        
        % ��֤���ʱ�̵�λ�����ٶ�ƫ�� (evaluate the position and velocity difference after injection)
        diffR = norm(Inject(2 : 4) - svTar(1 : 3));
        diffV = norm(Inject(5 : 7) - svTar(4 : 6));
        
        % ��飺���λ�ú��ٶ���� (check: whether the position and velocity difference meet the constraint)
        if diffR > aux.diffRMax || diffV > aux.diffVMax
            if aux.langFlag == 1
                disp(['��������',num2str(tarIndex),'��쾫�Ȳ�����!'])
            else
                disp(['Warning��for satellite ',num2str(tarIndex),', the injection error violates the constraint!'])
            end
            aux.errorFlag = 1;
        end
        
        % ��飺����ʣ������Լ�� (check: whether the mass after injection greater than dry mass)
        if Inject(8) < aux.mdry
            if aux.langFlag == 1
                disp(['��������',num2str(tarIndex),'��������С�ڸ���!'])
            else
                disp(['Warning��for satellite ',num2str(tarIndex),', the residual mass is lower than the dry mass!'])
            end
            aux.errorFlag = 1;
        end
        
        % ��ʾ��֤��� (display verification result of injection)
        if aux.errorFlag == 0
            if aux.langFlag == 1
                disp(['��ȷ����������', num2str(tarIndex),'�ɹ���죡'])
            else
                disp(['Correct��satellite ',num2str(tarIndex),' is injected into the orbit!'])
            end
        end
        
        % �ɹ�����visitAray��Ӧλ������ (if the satellite is injected into the target orbit successfully,
        % the corresponding element in the visitAray is replaced by 0)
        if aux.errorFlag == 0
            visitAray(tarIndex) = 0;
        else
            aux.errorFlag = 0;
        end
        
        % ���: ���ֲ���ʽ�ĵ����������Ƿ�����Լ�� (check: whether the number of satellites for two kind of deployments meet the constraint)
        switch typeFlag
            case 1 % �����»�� (for new launchings)
                if NsubTask > 16
                    if aux.langFlag == 1
                        disp('���󣺷����»����ʽ�ĵ�������������16��');
                    else
                        disp('Warning: for launching new rocket, the navigation satellites are greater than 16!');
                    end
                    aux.errorFlag = 1;
                end
            case 2 % ���ô��� (for piggyback missions)
                if NsubTask > 8
                    if aux.langFlag == 1
                        disp('���󣺴��ط�ʽ�ĵ�������������8��');
                    else
                        disp('Warning: for piggyback mission, the navigation satellites are greater than 8!');
                    end
                    aux.errorFlag = 1;
                end
        end
        
        % ��쵼�������� + 1 (number of injected satellites + 1)
        Nnavi = Nnavi + 1;
        
    end
    
    % ��ȡһ�� (read a new line)
    tline = fgets(fid);
    
end

% ���: ����������񱻶��ʹ�� (check: whether piggyback missions are used more than once)
if any(carryAray > 1)
    if aux.langFlag == 1
        disp('----------------------------------------------------------');
        disp(['���󣺴������񱻶��ʹ�ã����Ϊ: ' , num2str(find(carryAray > 1)')])
    else
        disp('----------------------------------------------------------');
        disp(['Warning��the piggyback missions are used more than once��indexes: ' , num2str(find(visitAray > 1)')])
    end
    aux.errorFlag = 1;
end

% ���: �Ƿ���δ������� (check: whether some satellites failded to inject into the target orbits)
if any(visitAray)
    if aux.langFlag == 1
        disp('----------------------------------------------------------');
        disp(['���󣺲�������δ����ȷ��죡���Ϊ: ' , num2str(find(visitAray ~= 0)')])
    else
        disp('----------------------------------------------------------');
        disp(['Warning��some satellites failed to inject into the target orbits��indexes: ' , num2str(find(visitAray ~= 0)')])
    end
    aux.errorFlag = 1;
end

% ����Obj2 (evaluate Obj2)
Obj2 = Nlaunch * aux.cLaunch + Ncarry * aux.cCarry + Nnavi * aux.cNavi;
if Obj2 > aux.minCost
    if aux.langFlag == 1
        disp(['������������ɱ�' , num2str(Obj2) , '����Լ��!']);
    else
        disp(['Warning��the cost of constellation violates the constraint!']);
    end
    aux.errorFlag = 1;
end

% ��ʾ'transfer.txt'����֤��� (display the verification result of 'transfer.txt')
if aux.errorFlag == 0
    if aux.langFlag == 1
        disp('transfer.txt ͨ����֤��')
    else
        disp('transfer.txt passed the verification��')
    end
end

% ������ (output verification results)
if aux.errorFlag == 0
    if aux.langFlag == 1
        disp('==========================================================');
        disp('ͨ��������֤�����������ƽ����')
        disp(['Obj1: ' , num2str(Obj1)])
        disp(['Obj2: ' , num2str(Obj2) , ' (���ҵ�λ)'])
        disp(['��������: ' , num2str(Nlaunch)])
        disp(['ʹ�ô�����: ' , num2str(Ncarry)])
        disp(['����������: ' , num2str(Nnavi)])
    else
        disp('==========================================================');
        disp('The results are listed below��')
        disp(['Obj1: ' , num2str(Obj1)])
        disp(['Obj2: ' , num2str(Obj2) , ' (currency unit)'])
        disp(['Number of new launches: ' , num2str(Nlaunch)])
        disp(['Number of piggybacks: ' , num2str(Ncarry)])
        disp(['Number of satellites: ' , num2str(Nnavi)])
    end
else
    if aux.langFlag == 1
        disp('==========================================================');
        disp('δ��ͨ��������֤��������ش���')
    else
        disp('==========================================================');
        disp('Failed to pass all validations��please check errors��')
    end
end

end


% =========================================================================
% �����������ṹ��
%
% �����
% aux              [structure]          ��������ṹ��
%
% ����: �ų�, �廪��ѧ�����
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Generate a structure of simulation parameters
%
% Output��
% aux              [structure]          structure of simulation parameters
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function aux = aux_
aux.mu = 398600; % ������������(km^3/s^2) (Gravitational constant of the Earth, km^3/s^2)
aux.req = 6378; % ����뾶(km) (Equatorial radius, km)
aux.j2 = 0.0010826; % j2�㶯(--) (Earth oblateness gravity coefficient, --)
aux.we = 7.29211585530e-5; % ������ת����(rad/sec) (Earth rotation rate, rad/sec)
aux.g0 = 9.80665; % ��ƽ���������ٶ�(m/s^2) (Gravity acceleration at sea level, m/s^2)
aux.Isp = 300; % �ȳ�(s) (Specific impulse, s)
aux.c = 1e-3 * aux.g0 * aux.Isp; % ȼ���ٶ� (Gas velocity, km/s)
aux.m0 = 100; % ���ǳ�ʼ����(kg) (Initial mass of navigation satellite, kg)
aux.mdry = 80; % ���Ǹ���(kg) (Dry mass of navigation satellite, kg)
aux.eleMin = 10 * pi / 180; % ������Ϳ�������(rad) (Minimum elevation angle, rad)
aux.dtr = pi / 180; % ��ת���� (deg to rad)
aux.rtd = 180 / pi; % ����ת�� (rad to deg)
aux.RelTol = 1e-13; % ��ֵ���־������ (Relative integration tolerance)
aux.AbsTol = 1e-13; % ��ֵ���������� (Absolute integration tolerance)
aux.diffTMax = 1e-3; % ��֤����ʱ�侫��(s) (time tolerance, s)
aux.diffRMax = 1e-3; % ��֤����λ�þ���(km) (position tolerance, km)
aux.diffVMax = 1e-6; % ��֤�����ٶȾ���(km/s) (velocity tolerance, km/s)
aux.diffMMax = 1e-4; % ��֤������������(kg) (mass tolerance, kg)
aux.tStep = 120; % ����׶β�������(s) (Sampling steps in service stage, s)
aux.JDtoMJD2000 = -2451544.5; % ��JD��MJD2000(day) (from JD to MJD2000, day)
aux.MJD2000toJD = 2451544.5; % ��MJD2000��JD(day) (from MJD2000 to JD, day)
aux.MJD2000_ti = jday(2020 , 1 , 1 , 0 , 0 , 0) + aux.JDtoMJD2000; % ����׶ο�ʼ����������(day) (epoch of construction begins, day)
aux.MJD2000_tf = jday(2020 , 4 , 1 , 0 , 0 , 0) + aux.JDtoMJD2000; % ����׶ν�������������(day) (epoch of construction ends, day)
aux.cLaunch = 1.2; % ���λ���ɱ�(���ҵ�λ) (Single launch cost, currency unit)
aux.cCarry = 0.2; % ���δ��سɱ�(���ҵ�λ) (Single piggyback cost, currency unit)
aux.cNavi = 0.05; % ���ŵ������ǳɱ�(���ҵ�λ) (Single navigation satellite cost, currency unit)
aux.minCost = 10; % ��������ɱ�(���ҵ�λ) (Constellation construction cost, currency unit)
aux.maxGDOP = 10; % ����GDOPԼ�� (GDOP constraint for each city)
aux.altiMin = 500; % ����߶�Լ�� (Minimum altitude constraint, km)
aux.altiPark = 900; % ������߶� (Altitude of parking orbit for launching new rockets, km)
aux.errorFlag = 0; % �����ʾ (Error flag during verification)
tt1 = [86400 * 0 : aux.tStep : 86400 * 1]; % ��1��(Sampling points of 1-th day, s)
tt2 = [86400 * 6 : aux.tStep : 86400 * 7]; % ��7��(Sampling points of 7-th day, s)
tt3 = [86400 * 29 : aux.tStep : 86400 * 30]; % ��30��(Sampling points of 30-th day, s)
aux.tAray = [tt1 , tt2 , tt3]; % ����׶ο�ʼ��ģ�����MJD2000=7396��ʼ������ʱ���, (Sampling points from MJD2000 = 7396, s)
end

% =========================================================================
% ����GDOP
%
% ����:
% satMtx          [nx3]          �ɼ�������ecfϵ�µ�λ������
% obsEcf          [1x3]          ������ecfϵ�µ�λ������
%
% ���:
% gdop            [1x1]           ���ξ�������(GDOP)
%
% ���ԣ�
% satMtx = [15524.471175 , -16649.826222 , 13512.2723887;
%          -2304.058534 , -23287.906465 , 11917.038105;
%           16680.243357 , -3069.625561 , 20378.551047;
%          -14799.931395 , -21425.35824 , 6069.947224];
% obsEcf = [-730.000 , -5440.000 , 3230.000];
% gdop_(satMtx , obsEcf)
%
% �ο����ף�
% https://en.wikipedia.org/wiki/Dilution_of_precision_(navigation)
%
% ����: �ų�, �廪��ѧ�����
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
% Example��
% satMtx = [15524.471175 , -16649.826222 , 13512.2723887;
%     -2304.058534 , -23287.906465 , 11917.038105;
%     16680.243357 , -3069.625561 , 20378.551047;
%     -14799.931395 , -21425.35824 , 6069.947224];
% obsEcf = [-730.000 , -5440.000 , 3230.000];
% gdop_(satMtx , obsEcf)
%
% Reference��
% https://en.wikipedia.org/wiki/Dilution_of_precision_(navigation)
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function gdop = getGDOP(satMtx , obsEcf)
seeNum = size(satMtx , 1); % �ɼ������� (number of visible satellites)
rTemp = satMtx - repmat(obsEcf , seeNum , 1); % ����վָ�����ǵ����� (vector from city to satellites)
rNorm = rTemp ./ repmat(sqrt(rTemp(: , 1).^2 + rTemp(: , 2).^2 + rTemp(: , 3).^2) , 1 , 3); % ����վָ�����ǵĵ�λ���� (unit vector from city to satellites)
H = [rNorm , ones(seeNum , 1)]; % ����H���� (generate H matrix)
gdop = sqrt( trace(pinv(H' * H)) ); % ����gdop (evaluate gdop)
end

% =========================================================================
% ����ָ��1
%
% ����:
% aux                  [structure]         �������
%    aux.MJD2000_tf    [1x1]               ����׶ν�������������(day)
%    aux.tAray         [1x2163]            ����׶ο�ʼ��ģ�����MJD2000=7396��ʼ������ʱ���(s)
%    aux.city          [structure]         ������Ϣ�ṹ��
%    aux.cons          [structure]         ������Ϣ�ṹ��
%    aux.eleMin        [1x1]               ������Ϳ�������(rad)
%    aux.maxGDOP       [1x1]               ����GDOPԼ��
%
% ���:
% Obj1                 [1x1]               ���ָ��1
%
% ����: �ų�, �廪��ѧ�����
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
Nt = size(aux.tAray , 2); % ʱ��������� (number of sampling points)
Ncity = size(aux.city.ecf , 1); % ������ (number of cities)
NnaviSat = size(aux.cons.oe , 1); % ���������� (number of navigation satellites)
eleMin = aux.eleMin; % ������Ϳɼ����� (minimum elevation angle)

% ����������ecfϵ�µ�λ���ٶ� (evaluate position and velocity under ECF frame)
consEcfMtx = zeros(NnaviSat * 6 , Nt);
count = 1;
for i = 1 : NnaviSat % ����ѭ�� (loop of navigation satellites)
    % ����ecf��λ���ٶ� (evaluate position and velocity under ECF frame)
    [rv_eci , rv_ecf] = oe2eciecf_(aux.cons.oe(i , :) , aux.MJD2000_tf , aux.MJD2000_tf + aux.tAray / 86400 , aux);
    % ��������λ���ٶ� (save position and velocity under ECF frame)
    consEcfMtx( (count - 1) * 6 + 1 : count * 6 , : ) = rv_ecf;
    count = count + 1;
end

% �������г��������в������ϵ�GDOPֵ (evaluate GDOP matrix for each city and each sampling point)
GDOPMtx = 1e8 * ones(Ncity , Nt); % ��ʼ��GDOP���� (pre allocate GDOP matrix)

% ʱ��ѭ�� (time loop)
for j = 1 : Nt
    
    % ����ѭ�� (city loop)
    for k = 1 : Ncity
        
        % ��ȡ��k����������, ECFϵ (extract position of the k-th city, ECF frame)
        cityEcf = aux.city.ecf(k , :);
        
        % ��ȡ�������������ڵ�j���������λ���ٶȣ�ECFϵ (extract position of constellation at the j-th sampling point, ECF frame)
        satEcf = reshape(consEcfMtx(: , j) , 6 , NnaviSat)';
        
        % ����ӵ�k�����е������������ǵ�λ��ʸ�� (evaluate vectors from the k-th city to all satellites of constellation)
        city2sat = [satEcf(: , 1) - cityEcf(1) , satEcf(: , 2) - cityEcf(2) , satEcf(: , 3) - cityEcf(3)];
        
        % ����������Գ��еļн� (evalute angles between cityEcf vector and city2sat vectors)
        thetaAray = acos( (cityEcf(1) .* city2sat(: , 1) + cityEcf(2) .* city2sat(: , 2) + cityEcf(3) .* city2sat(: , 3)) ./ ...
            ( sqrt(cityEcf(1)^2 + cityEcf(2)^2 + cityEcf(3)^2) * sqrt(city2sat(: , 1).^2 + city2sat(: , 2).^2 + city2sat(: , 3).^2) ) );
        
        % ��ȡ�ɼ����Ǳ�� (extract visible satellite indexes)
        seeIndex = find(thetaAray < (0.5 * pi - eleMin));
        
        % ����ɼ������� >= 4������GDOP (if the number of visiable satellties is greater than 4, evaluate GDOP)
        if size(seeIndex , 1) >= 4
            GDOPMtx(k , j) = getGDOP(satEcf(seeIndex , 1 : 3) , cityEcf);
        end
        
    end
    
end

% ����Ŀ�꺯�� (evaluate Obj1)
Obj1 = ( max(GDOPMtx') <= aux.maxGDOP ) * aux.city.weight;

end


% =========================================================================
% ͨ��ƫ���ʺ�ƽ����Ǽ���ƫ����Ǻ�������
%
% ����:
% e            [1x1]         ƫ����
% M            [nx1]         ƽ�����
%
% ���:
% Ek1          [nx1]         ƫ�����
% theta        [nx1]         ������
%
% ����: �ų�, �廪��ѧ�����
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
% J2�㶯�������ѧģ��
%
% ����:
% t             [1x1]            ʱ��
% x             [6x1]            ״̬
% aux           [structure]      �������
%   aux.mu      [1x1]            ������������(km^3/s^2)
%   aux.req     [1x1]            ����뾶(km)
%   aux.j2      [1x1]            j2�㶯(--)
%
% ���:
% dx            [6x1]            ������
%
% ����: �ų�, �廪��ѧ�����
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
% Dynamic model under j2 perturbation
%
% Input��
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
% �������תeci��ecf�µ�λ���ٶ�
%
% ����:
% xx_oe             [1x6]           t0ʱ�̵Ĺ������
%                                   [sma(km) , ecc(--) , incl(rad) , ...
%                                    raan(rad) , argper(rad) , manom(rad)]
% t0                [1x1]           xx_oe����Ԫ, MJD2000 (day)
% tt                [1xn]           ��Ҫ�����ʱ���, MJD2000 (day)
% aux               [structure]     �������
%   aux.mu          [1x1]           ������������(km^3/s^2)
%   aux.req         [1x1]           ����뾶(km)
%   aux.j2          [1x1]           j2�㶯(--)
%   aux.we          [1x1]           ������ת����(rad/sec)
%   aux.MJD2000toJD [1x1]           ��MJD2000��JD(day)
%
% �����
% rv_eci            [6xn]           ttʱ��eciλ���ٶ� [x,y,z,vx,vy,vz] (km, km/s)
% rv_ecf            [6xn]           ttʱ��ecfλ���ٶ� [x,y,z,vx,vy,vz] (km, km/s)
%
% ����: �ų�, �廪��ѧ�����
% chenzhang86@tsinghua.edu.cn
% 2017/09/01 (Sep/01/2017)
% -------------------------------------------------------------------------
% Orbital elements to position and velocity under both ECI and ECF frame
%
% Input��
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
% Output��
% rv_eci            [6xn]           states under ECI frame at sampling points tt [x,y,z,vx,vy,vz] (km, km/s)
% rv_ecf            [6xn]           states under ECF frame at sampling points tt [x,y,z,vx,vy,vz] (km, km/s)
%
% Author: Chen ZHANG, Tsinghua Space Center
% chenzhang86@tsinghua.edu.cn
% Sep/01/2017
% -------------------------------------------------------------------------
function [rv_eci , rv_ecf] = oe2eciecf_(xx_oe , t0 , tt , aux)

% ��ȡ���� (extract simulation parameters)
mu = aux.mu;
j2 = aux.j2;
req = aux.req;
we = aux.we;

% ��ȡt0ʱ�̵Ĺ������ (extract orbital elements at t0)
sma = xx_oe(1); % �볤�� (km) (semi-major axis, km)
ecc = xx_oe(2); % ƫ���� (--) (eccentricity, --)
incl = xx_oe(3); % ������ (rad) (inclination, rad)
raan0 = xx_oe(4); % ������ྭ(rad) (right ascension of ascending node, rad)
argper0 = xx_oe(5); % ���ص���� (rad) (argument of perigee, rad)
manom0 = xx_oe(6); % ƽ����� (rad) (mean anomaly, rad)

tt = tt(:); % ��Ҫ�����ʱ��tt��Ϊ������ (transform tt to column vector)

% ����j2������ƽ���仯�� (evaluate the change rate of elements under j2 perturbation)
t1 = 1 - ecc.^2;
t2 = (req ./ sma).^2;
An = sqrt(mu ./ sma.^3);
draan = -1.5 * An * j2 ./ t1.^2 .* t2 .* cos(incl);
dargper = 1.5 * An * j2 ./ t1.^2 .* t2 .* (2 - 2.5 * sin(incl).^2);
dmanom = 0.75 * An * j2 ./ t1.^(1.5) .* t2 .* (3 * cos(incl).^2 - 1);

% ��������㴦��ƽ�� (evaluate the orbital elements at sampling points)
raan = mod(raan0 + draan .* (tt - t0) * 86400 , 2 * pi);
argper = mod(argper0 + dargper .* (tt - t0)  * 86400 , 2 * pi);
manom = mod(manom0 + (An + dmanom) .* (tt - t0)  * 86400 , 2 * pi);

% ���������� (evaluate true anomaly)
[~ , tanom] = m2ef_(ecc , manom);

% ����eciλ���ٶ� (evaluate position and velocity under ECI frame)
r = sma .* (1 - ecc .^ 2) ./ (1 + ecc .* cos(tanom));
v = sqrt( 2 * mu ./ r - mu ./ sma );
gamma = atan( ecc .* sin(tanom) ./ (1 + ecc .* cos(tanom)) );
x = r .* ( cos(tanom + argper) .* cos(raan) - sin(tanom + argper) .* cos(incl) .* sin(raan) );
y = r .* ( cos(tanom + argper) .* sin(raan) + sin(tanom + argper) .* cos(incl) .* cos(raan) );
z = r .* ( sin(tanom + argper) .* sin(incl) );
vx = v .* ( -sin(tanom + argper - gamma) .* cos(raan) - cos(tanom + argper - gamma) .* cos(incl) .* sin(raan) );
vy = v .* ( -sin(tanom + argper - gamma) .* sin(raan) + cos(tanom + argper - gamma) .* cos(incl) .* cos(raan) );
vz = v .* ( cos(tanom + argper - gamma) .* sin(incl));

% ����eci�µ�λ���ٶ� (save position and velocity under ECI frame)
r_eci = [x , y , z]';
v_eci = [vx , vy , vz]';
rv_eci = [r_eci ; v_eci];

% eciתecf (transform eci to ecf)
if nargout > 1
    % ���������tt�ĸ������κ���ʱ(rad) (evaluate Greenwich sidereal time, rad)
    vgast = gstime(tt + aux.MJD2000toJD);
    
    % ����������ת���� (evaluate coordinate transformation matrix)
    tm = zeros(3 , 3 , length(tt));
    tmdot = zeros(3 , 3 , length(tt));
    tm(1 , 1 , :) = cos(vgast); tm(1 , 2 , :) = sin(vgast); tm(1 , 3 , :) = 0;
    tm(2 , 1 , :) = -sin(vgast); tm(2 , 2 , :) = cos(vgast); tm(2 , 3 , :) = 0;
    tm(3 , 1 , :) = 0; tm(3 , 2 , :) = 0; tm(3 , 3 , :) = 1;
    tmdot(1 , 1 , :) = -we * sin(vgast); tmdot(1 , 2 , :) = we * cos(vgast); tmdot(1 , 3 , :) = 0;
    tmdot(2 , 1 , :) = -we * cos(vgast); tmdot(2 , 2 , :) = -we * sin(vgast); tmdot(2 , 3 , :) = 0;
    tmdot(3 , 1 , :) = 0; tmdot(3 , 2 , :) = 0; tmdot(3 , 3 , :) = 0;
    
    % ����ecf�µ�λ���ٶ� (evaluate position and velocity under ECF frame)
    r_ecf = x_multi_y_(tm , r_eci);
    v_ecf = x_multi_y_(tm , v_eci) + x_multi_y_(tmdot , r_eci);
    rv_ecf = [r_ecf ; v_ecf];
end

end

% =========================================================================
% J2ƽ����������������epochAray��λ���ٶȺ͹������
%
% ����:
% sat                 [structure]         ���ǲ����ṹ��
%    sat.id           [nx1]               ���Ǳ��
%    sat.MJD2000      [nx1]               ������Ԫ
%    sat.oe           [nx6]               ���ǹ������
% epochAray           [1xm]               ��Ҫ�����ʱ��
% aux                 [structure]         �������
%    aux.mu           [1x1]               ������������(km^3/s^2)
%    aux.req          [1x1]               ����뾶(km)
%    aux.j2           [1x1]               j2�㶯(--)
%
% ���:
% svCubic             [nx6xm]             ����������epochAray��λ���ٶ�
% oeCubic             [nx6xm]             ����������epochAray�Ĺ������
%
% ����: �ų�, �廪��ѧ�����
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

% ��ȡ������� (extract simulation parameters)
mu = aux.mu;
j2 = aux.j2;
req = aux.req;

% ��ȡ��Ҫ�����ʱ�������������� (extract number of sampling points and number of satellites)
epochNum = size(epochAray , 2);
satNum = size(sat.oe , 1);

% ��ȡ�������ǵĲ��� (extract epoch and orbital elements of all satellites)
epoch = sat.MJD2000;
sma = sat.oe(: , 1);
ecc = sat.oe(: , 2);
incl = sat.oe(: , 3);
raan0 = sat.oe(: , 4);
argper0 = sat.oe(: , 5);
manom0 = sat.oe(: , 6);

% ����j2������ƽ���仯�� (evaluate the change rate of elements under j2 perturbation)
t1 = 1 - ecc.^2;
t2 = (req ./ sma).^2;
An = sqrt(mu ./ sma.^3);
draan = -1.5 * An * j2 ./ t1.^2 .* t2 .* cos(incl);
dargper = 1.5 * An * j2 ./ t1.^2 .* t2 .* (2 - 2.5 * sin(incl).^2);
dmanom = 0.75 * An * j2 ./ t1.^(1.5) .* t2 .* (3 * cos(incl).^2 - 1);

svCubic = zeros(satNum , 6 , epochNum); % λ���ٶ�Ԥ���� (pre allocate space for position and velocity)
oeCubic = zeros(satNum , 6 , epochNum); % �������Ԥ���� (pre allocate space for orbital elements)
for i = 1 : epochNum
    
    % ��ȡ��Ҫ�����ʱ�� (extract sampling points)
    epochNow = epochAray(i);
    
    % ��������㴦��ƽ�� (evaluate orbital elements at sampling points)
    raan = mod(raan0 + draan .* (epochNow - epoch) * 86400 , 2 * pi);
    argper = mod(argper0 + dargper .* (epochNow - epoch) * 86400 , 2 * pi);
    manom = mod(manom0 + (An + dmanom) .* (epochNow - epoch) * 86400 , 2 * pi);
    
    % ���������� (evaluate true anomaly)
    [~ , tanom] = m2ef_(ecc , manom);
    
    % ����λ���ٶ� (evaluate position and velocity)
    r = sma .* (1 - ecc.^2) ./ (1 + ecc .* cos(tanom));
    v = sqrt(2 * mu ./ r - mu ./ sma );
    gamma = atan(ecc .* sin(tanom) ./ (1 + ecc .* cos(tanom)));
    x = r .* (cos(tanom + argper) .* cos(raan) - sin(tanom + argper) .* cos(incl) .* sin(raan));
    y = r .* (cos(tanom + argper) .* sin(raan) + sin(tanom + argper) .* cos(incl) .* cos(raan));
    z = r .* (sin(tanom + argper) .* sin(incl) );
    vx = v .* (-sin(tanom + argper - gamma) .* cos(raan) - cos(tanom + argper - gamma) .* cos(incl) .* sin(raan));
    vy = v .* (-sin(tanom + argper - gamma) .* sin(raan) + cos(tanom + argper - gamma) .* cos(incl) .* cos(raan));
    vz = v .* (cos(tanom + argper - gamma) .* sin(incl));
    
    % ����λ���ٶȣ�������� (save position velocity and orbital elements)
    svCubic(: , : , i) = [x , y , z , vx , vy , vz];
    oeCubic(: , : , i) = [sma , ecc , incl , raan , argper , manom];
end

end

% =========================================================================
% J2ƽ������sat�ṹ���б��Ϊindex�����ǣ��õ�Ŀ����Ԫʱ��(tarEpoch)��λ��
% �ٶȺ͹������
%
% ����:
% sat                  [structure]        ���ǲ���
%    sat.id            [nx1]              ���Ǳ��
%    sat.MJD2000       [nx1]              ���ǳ�ʼ��Ԫʱ��
%    sat.oe            [nx6]              �����ڳ�ʼ��Ԫ�Ĺ������
%    sat.sv            [nx6]              �����ڳ�ʼ��Ԫ��λ���ٶ�
% tarEpoch             [1x1]              ��Ҫ�������Ԫʱ��
% aux                  [structure]        �������
%    aux.mu            [1x1]              ������������(km^3/s^2)
%    aux.req           [1x1]              ����뾶(km)
%    aux.j2            [1x1]              j2�㶯(--)
% index                [1x1]              ��Ҫ��������Ǳ��
%
% ���:
% sv                   [1x6]              ������tarEpoch��λ���ٶ�
% oe                   [1x6]              ������tarEpoch�Ĺ������
%
% ����: �ų�, �廪��ѧ�����
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

% ��ȡ������� (extract simulation parameters)
mu = aux.mu;
j2 = aux.j2;
req = aux.req;

% ��ȡsat�б��Ϊindex�������ڳ�ʼ��Ԫʱ�̵Ĺ������
epoch = sat.MJD2000(index);
sma = sat.oe(index , 1);
ecc = sat.oe(index , 2);
incl = sat.oe(index , 3);
raan0 = sat.oe(index , 4);
argper0 = sat.oe(index , 5);
manom0 = sat.oe(index , 6);

% ��ȡ��Ҫ�������Ԫʱ�� (extract the target epoch)
tarEpoch = tarEpoch(:);
epochNum = size(tarEpoch , 1);

% ����j2������ƽ���仯�� (evaluate the change rate of elements under j2 perturbation)
t1 = 1 - ecc.^2;
t2 = (req ./ sma).^2;
An = sqrt(mu ./ sma.^3);
draan = -1.5 * An * j2 ./ t1.^2 .* t2 .* cos(incl);
dargper = 1.5 * An * j2 ./ t1.^2 .* t2 .* (2 - 2.5 * sin(incl).^2);
dmanom = 0.75 * An * j2 ./ t1.^(1.5) .* t2 .* (3 * cos(incl).^2 - 1);

% ����tarEpoch����ƽ�� (evaluate orbital elements at sampling points)
raan = mod(raan0 + draan .* (tarEpoch - epoch) * 86400 , 2 * pi);
argper = mod(argper0 + dargper .* (tarEpoch - epoch) * 86400 , 2 * pi);
manom = mod(manom0 + (dmanom + An) .* (tarEpoch - epoch) * 86400 , 2 * pi);

% ���������� (evaluate true anomaly)
[~ , theta] = m2ef_(ecc , manom);

% ����λ���ٶ� (evaluate position and velocity)
r = sma .* (1 - ecc.^2) ./ (1 + ecc .* cos(theta));
v = sqrt(2 * mu ./ r - mu ./ sma );
gamma = atan(ecc .* sin(theta) ./ (1 + ecc .* cos(theta)));
x = r .* (cos(theta + argper) .* cos(raan) - sin(theta + argper) .* cos(incl) .* sin(raan));
y = r .* (cos(theta + argper) .* sin(raan) + sin(theta + argper) .* cos(incl) .* cos(raan));
z = r .* (sin(theta + argper) .* sin(incl) );
vx = v .* (-sin(theta + argper - gamma) .* cos(raan) - cos(theta + argper - gamma) .* cos(incl) .* sin(raan));
vy = v .* (-sin(theta + argper - gamma) .* sin(raan) + cos(theta + argper - gamma) .* cos(incl) .* cos(raan));
vz = v .* (cos(theta + argper - gamma) .* sin(incl));

% ����λ���ٶȣ�������� (save position velocity and orbital elements)
sv = [x , y , z , vx , vy , vz];
oe = [sma * ones(epochNum , 1) , ecc * ones(epochNum , 1), incl * ones(epochNum , 1) , raan , argper , manom];

end

% -------------------------------------------------------------------------
% ��黬�жεĻ��־���
%
% ����:
% state0                    [1x8]               ���п�ʼ״̬
% state1                    [1x8]               ���н���״̬
% aux                       [structure]         �������
%    aux.RelTol             [1x1]               ��ֵ���־������
%    aux.RelTol             [1x1]               ��ֵ����������
%
% ���:
% diffR                     [1x1]               λ�����
% diffV                     [1x1]               �ٶ����
%
% ����: �ų�, �廪��ѧ�����
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
% ��֤����Ƿ�����
%
% ����:
% coast0                    [1x8]                 ��һ��״̬
% xcraft                    [1x8]                 �ڶ���״̬
% aux                       [structure]           �������
%   aux.diffTMax            [1x1]                 ��֤����ʱ�侫�� (s)
%   aux.diffRMax            [1x1]                 ��֤����λ�þ��� (km)
%   aux.diffVmax            [1x1]                 ��֤�����ٶȾ��� (km/s)
%   aux.diffMMax            [1x1]                 ��֤������������ (kg)
%   aux.langFlag            [1x1]                 ���Ա�ʶ
%   aux.errorFlag           [1x1]                 �����ʶ
%
% ���:
% aux                       [structure]           �������
%
% Copyright �ų�(Chen Zhang)
% �廪��ѧ�����(Tsinghua Space Center)
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
        disp('���󣺵������ǵ�״̬��������')
    else
        disp(['Warning: the states are not continuous!']);
    end
    aux.errorFlag = 1;
end
end

% =========================================================================
% �þ��������������ת��
%
% ����:
% x           [3x3xn]       ����ת������
% y           [3xn]         λ�û��ٶȾ���
%
% ���:
% x_multi_y   [3xn]         ����ת�����
%
% ����: �ų�, �廪��ѧ�����
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
