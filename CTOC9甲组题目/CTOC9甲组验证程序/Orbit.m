function [TH, RVH, Tobv0, Tobvf, RVobv0] = Orbit(T0, RV0, m0, Sat, options, TaskNo, SatNo)
%==========================================================================
%   ���Ө    yangkeying@bit.edu.cn      2017-09-10
%==========================================================================
% ���ܣ�ʵ�ֶַι�����֣����λ�á��ٶȡ��������ȣ���ȡ����ͼ��״̬����Ϣ
% ���룺 
%      T0           �������ʼʱ�̣���λ��s
%      RV0          ������ǵĳ�ʼ״̬��6ά����������λ��km, km/s
%      m0           ������ǵĳ�ʼ��������λ��kg
%      Sat          �������иü�����ǵ����ݾ���N*17ά����
%      options      ����������
%      TaskNo       ��ǰ���������
%      SatNo        �ü���������
% �����
%      TH           ��ǰ�������иü�����ǵ�ʱ������㣬��λ��s
%      RVH          ��ǰ�������иü������״̬��ʱ�����̣���λ��km, km/s
%      Tobv0        �ü�����ǽ������ʱ�̣���λ��s
%      Tobvf        �ü�����ǽ�������ʱ�̣���λ��s
%      RVobv0       �ü�����ǽ�����ʱ��״̬��6ά����������λ��km, km/s
%==========================================================================

global Re mu Isp ge

TH = [];
RVH = [];
Count1 = 0; % ���ڼ�¼״̬1���ֵĴ���
Count2 = 0; % ���ڼ�¼״̬2���ֵĴ���
Count3 = 0; % ���ڼ�¼״̬3���ֵĴ���

% ���ʱ�䵥����
Ttemp = Sat(:,3);
if any(diff(Ttemp)<0)
    fprintf('�� %d ���������ʱ�䲻�ǵ�������\n',TaskNo);
end

% ���л���
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
        fprintf('event��ȡֵֻ���� 1, 2, 3, 4.\n');
    end
    
    % ����ǰ����ε�ĩ״̬��ֵΪ��һ����εĳ�״̬
    T0 = Sat(jj,3);
    RV0 = Sat(jj,5:10);
    m0 = Sat(jj,11);
    
    % ��¼��������������ʱ�̺�״̬ʱ������
    TH = cat(1, TH, tt(1:end-1));
    RVH = cat(1, RVH, yy(1:end-1,1:6));
end
if Count1 > 1 || Sat(end-1,15) ~= 1
    fprintf('�� %d ���������� ���� %d ��״̬1���ֵĴ�����λ������\n',TaskNo, SatNo);
end
if Count2 > 1 || Sat(end,15) ~= 2
    fprintf('�� %d ���������� ���� %d ��״̬2���ֵĴ�����λ������\n',TaskNo, SatNo);
end
if Count3 > 5
    fprintf('�� %d ���������� ���� %d �Ļ������������� 5 ��\n',TaskNo, SatNo);
end


function testRVM(TaskNo, SatNo, ColNo, rf, rt, vf, vt, mf, mt)
% ���ڼ��λ�á��ٶȡ���������

global rTol vTol mTol

if abs(norm(rf)-norm(rt)) > rTol
    fprintf('�ڵ� %d ���������У��� %d �ż�������ڵ� %d �е�λ�þ��Ȳ�����Ҫ��,λ�����Ϊ %8f (km)\n',TaskNo, SatNo, ColNo,abs(norm(rf)-norm(rt)))
end
if abs(norm(vf)-norm(vt)) > vTol
    fprintf('�ڵ� %d ���������У��� %d �ż�������ڵ� %d �е��ٶȾ��Ȳ�����Ҫ��,�ٶ����Ϊ %8f (km/s)\n',TaskNo, SatNo, ColNo,abs(norm(vf)-norm(vt)))
end
if abs(norm(mf)-norm(mt)) > mTol
    fprintf('�ڵ� %d ���������У��� %d �ż�������ڵ� %d �е��������Ȳ�����Ҫ��,�������Ϊ %8f (km)\n',TaskNo, SatNo, ColNo,abs(norm(mf)-norm(mt)))
end