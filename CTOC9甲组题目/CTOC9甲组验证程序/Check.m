%==========================================================================
%               �ھŽ�ȫ���ռ�����ƾ����������������
%==========================================================================
%   ����     qirui@bit.edu.cn      2017-09-10
%==========================================================================

%��չ�����
clc;clear

%==========================================================================
%���Ȳ�����Ҫ�󣬸�������ж�Ϊ���Υ��Լ�����ض�Υ��Լ�����ж�ԭ�������ⷽ����
%����2���ʱ������޸ģ��޸ĺ�÷ֲ��ܳ���ԭ�����÷֣�������������ԭ��������һ��
%==========================================================================

%% =========================��������===============================
%��ʼ��ȫ�ֱ���
global mu J2 Re ge Isp 
global rTol vTol mTol

mu   = 398600;    %������������ km2/s2
J2   = 0.0010826; %J2�㶯����
Re   = 6378.137;  %����뾶 km
ge   = 9.80665*1e-3;   %����ƽ����������ٶȣ� km/sec^2
Isp  = 300;       %�ƽ����ȳ壬��λ��s
SPstep = 1;       %����Ǽ�����ʱ�������������λ��s

% λ�á��ٶȡ��������������㾫��
rTol = 1e-3; % ��λ��km
vTol = 1e-6; % ��λ��km/s
mTol = 1e-4; % ��λ��kg
GainTol = 1e-3; % ��λ��km^2

% �������������Ͳ���
RelTol = 1e-13 ; 
AbsTol = 1e-13;
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

GEOgain = zeros(18,1);

%% ==========================���ݶ���==============================
% ����������
data = load('result_test.txt');

% �ٷ�GEO���Ǿ��ȡ����������ǡ���λ��
GEO = load('GEO.txt');


%% ========================��ʼ״̬���============================
% ��ʼ���Ҫ�غ�����
sma0 = 6778.137; %�볤�ᣬ��λ��km
ecc0 = 0; %ƫ����
inc0 = 42.8*pi/180; %�����ǣ���λ��rad
m0 = 500; %���������������λ��kg

% �������ż�����ǵĳ�ʼ�������
fprintf('================================��ʼ״̬���==================================== \n\n')
for kk = 1:3
    flagIni = 0;
    rv0(kk,:)  = data(min(find(data(:,4) == kk)),5:10);
    ele0(kk,:) = RV2Ele(rv0(kk,:),mu);
    dsma = ele0(kk,1) - sma0;
    decc = ele0(kk,2) - ecc0;
    dinc = ele0(kk,3) - inc0;
    dm   = data(min(find(data(:,4) == kk)),11) - m0;
    if abs(dsma) > 1e-5
        fprintf('��%d�ż������Υ����ʼ�볤��Լ��������볤��%f��ʵ�ʰ볤��%f\n\n',kk, sma0, ele0(kk,1))
        flagIni = flagIni + 1;
    end
    if abs(decc) > 1e-7
        fprintf('��%d�ż������Υ����ʼƫ����Լ��������ƫ����%f��ʵ��ƫ����%f\n\n',kk, ecc0, ele0(kk,2))
        flagIni = flagIni + 1;
    end
    if abs(dinc) > 1e-7
        fprintf('��%d�ż������Υ����ʼ������Լ�����������%f��ʵ�����%f\n\n',kk, inc0, ele0(kk,3))
        flagIni = flagIni + 1;
    end
    if abs(dm) > 1e-10
        fprintf('��%d�ż������Υ����ʼ����Լ������������%f��ʵ������%f\n\n',kk, m0, data(min(find(data(:,4) == kk)),11))
        flagIni = flagIni + 1;
    end
    if flagIni == 0
        fprintf('��%d�ż�����ǳ�ʼ״̬��֤�ɹ� \n\n', kk)
    else
        fprintf('��%d�ż�����ǳ�ʼ״̬��֤ʧ�� \n\n', kk)
    end
end
fprintf('=============================================================================== \n\n')

%% =============������֣����ȼ�⣬��ȡ�����Ϣ����������=============
% ��ʼ״̬
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

for ii = 1 : data(end,2) % data(end,2)�����������
    
    % ���ݶ�ȡ
    SubTask =  data(find(data(:,2) == ii),:); % ��ȡ������ii���������ݣ�����SubTask
    Sat1 = SubTask(find(SubTask(:,4) == 1),:); % ��ȡ������ii�и�������ǵ��������ݣ��ֱ𸳸�Sat1��Sat2��Sat3
    Sat2 = SubTask(find(SubTask(:,4) == 2),:);
    Sat3 = SubTask(find(SubTask(:,4) == 3),:);
    
    if Sat1(find(Sat1(:,15)==1),3) == T1End
        disp('�������1����һ�������������ʱ�䲻���뵱ǰ�����������ʱ����ͬ')
    end
    if Sat2(find(Sat2(:,15)==1),3) == T2End
        disp('�������2����һ�������������ʱ�䲻���뵱ǰ�����������ʱ����ͬ')
    end
    if Sat3(find(Sat3(:,15)==1),3) == T3End
        disp('�������3����һ�������������ʱ�䲻���뵱ǰ�����������ʱ����ͬ')
    end
    
    %Orbit����ʵ�ֶַι�����֣����λ�á��ٶȡ��������Ⱥ͹���߶ȣ���ȡ����ͼ��״̬����Ϣ
    [TH1, RVH1, Tobv10, Tobv1f, RVobv10] = Orbit(T10, RV10, m10, Sat1, options, ii, 1);
    if any((sqrt(RVH1(:,1).^2 + RVH1(:,2).^2 + RVH1(:,3).^2) - Re) <300) == 1
        fprintf('�� %d ���������� ���� %d ����߶�С�� 300 km\n', ii, 1);
    end
    
    [TH2, RVH2, Tobv20, Tobv2f, RVobv20] = Orbit(T20, RV20, m20, Sat2, options, ii, 2);
    if any((sqrt(RVH2(:,1).^2 + RVH2(:,2).^2 + RVH2(:,3).^2) - Re) <300) == 1
        fprintf('�� %d ���������� ���� %d ����߶�С�� 300 km\n', ii, 2);
    end
    
    [TH3, RVH3, Tobv30, Tobv3f, RVobv30] = Orbit(T30, RV30, m30, Sat3, options, ii, 3);
    if any((sqrt(RVH3(:,1).^2 + RVH3(:,2).^2 + RVH3(:,3).^2) - Re) <300) == 1
        fprintf('�� %d ���������� ���� %d ����߶�С�� 300 km\n', ii, 1);
    end
    
    % �����������ǽ�������ʱ���һ����
    if Tobv10 ~= Tobv20 || Tobv10 ~= Tobv30
        fprintf('�� %d �������������ż�����ǽ�������ʱ�䲻һ��\n', ii);
    end
    if Tobv1f ~= Tobv2f || Tobv1f ~= Tobv3f
        fprintf('�� %d �������������ż�����ǳ�������ʱ�䲻һ��\n', ii);
    end
    
    % ������ż�����Ǽ����Ծ����Ƿ�����Լ��
    Tsp = TH1(1) : SPstep : TH1(end);
    RVsp1 = spline(TH1', RVH1', Tsp);
    RVsp2 = spline(TH2', RVH2', Tsp);
    RVsp3 = spline(TH3', RVH3', Tsp);
    for gg = 1 : length(Tsp)
        if norm(RVsp1(1:3, gg) - RVsp2(1:3, gg)) < 0.1
            fprintf('�� %d ���������е� 1 �� �� 2 �ż������֮��ľ����� t = %f ʱС��100m\n', ii, Tsp(gg));
        elseif norm(RVsp2(1:3, gg) - RVsp3(1:3, gg)) < 0.1
            fprintf('�� %d ���������е� 2 �� �� 3 �ż������֮��ľ����� t = %f ʱС��100m\n', ii, Tsp(gg));
        elseif norm(RVsp1(1:3, gg) - RVsp3(1:3, gg)) < 0.1
            fprintf('�� %d ���������е� 1 �� �� 3 �ż������֮��ľ����� t = %f ʱС��100m\n', ii, Tsp(gg));
        end
    end
 
    % ��ȡ�����ڸ��������������ʱ��״̬     
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
    
    % �������棬����������Ǿ��ڲ����ڣ����ڵ���ͬ��
    GEONo  = Sat1(end,16);
    lamda  = GEO(GEONo,2);
    El     = GEO(GEONo,3);
    Az     = GEO(GEONo,4);    
    GainComp = Gain(ii, ttN1, yyN1, yyN2, yyN3, lamda, El, Az); % ������ĵ�ǰ�������õ�GEO�������  
    GainProv = (SubTask(end,17) - SubTask(1,17)); % �������ṩ�ĵ�ǰ�������õ�GEO�������  
    if abs(GainComp - GainProv) > GainTol
        fprintf('�� %d ��������ļ�����治��ȷ,����������Ϊ %8f (km^2)\n',ii,abs(GainComp - GainProv));
    end
    GEOgain(GEONo) = GEOgain(GEONo) + GainProv;
    
    % ��ǰ�������ĩ״̬��ֵΪ��һ��������ĳ�״̬

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


    fprintf('�� %d / %d ����������֤���� \n\n', ii, data(end,2)) % data(end,2)�����������
end


%% =========================������===============================
if min(GEOgain) >= 1e7
    Obj1 = Tobv1f;
    mf = [Sat1(end,11);Sat2(end,11);Sat3(end,11)];
    if min(m0-mf) > 200
        error('ȼ�����ĳ��� 200 kg');
    end
    Obj2 = sum(m0 - mf);
    disp('����������ɣ�')
    
    % ������ָ��
    disp('��һָ�� (day)')
    disp(Obj1/3600/24) % �����һָ�꣬��λ��day
    disp('�ڶ�ָ�� (kg)')
    disp(Obj2)
    disp('���������ȼ������ (kg)')
    disp(m0 - mf)
        
    % �����GEO���ǵ�����
    disp('GEO��������')
    disp(GEOgain)

else
    disp('������δ��ɣ�')
end