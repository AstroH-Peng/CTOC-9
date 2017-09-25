% =========================================================================
% 计算GDOP
%
% 输入:
% satMtx    [n x 3]    可见卫星在ecf系下的位置坐标
% obsEcf    [1 x 3]    城市在ecf系下的位置坐标
%
% 输出:
% gdop      [1 x 1]     几何精度因子
%
% 测试：
% satMtx = [15524.471175 , -16649.826222 , 13512.2723887;
%     -2304.058534 , -23287.906465 , 11917.038105;
%     16680.243357 , -3069.625561 , 20378.551047;
%     -14799.931395 , -21425.35824 , 6069.947224];
% obsEcf = [-730.000 , -5440.000 , 3230.000];
% gdop_(satMtx , obsEcf)
%
% 参考文献：
% https://en.wikipedia.org/wiki/Dilution_of_precision_(navigation)
%
% Copyright 张晨(Chen Zhang)
% 清华大学宇航中心(Tsinghua Space Center)
% chenzhang86@tsinghua.edu.cn
% 2017/09/01
% -------------------------------------------------------------------------
function gdop = getGDOP(satMtx , obsEcf)
seeNum = size(satMtx , 1); % 可见卫星数
rTemp = satMtx - repmat(obsEcf , seeNum , 1); % 地面站指向卫星的向量
rNorm = rTemp ./ repmat(sqrt(rTemp(: , 1).^2 + rTemp(: , 2).^2 + rTemp(: , 3).^2) , 1 , 3); % 地面站指向卫星的单位向量
H = [rNorm , ones(seeNum , 1)]; % 构造H矩阵
gdop = sqrt( trace(pinv(H' * H)) ); % 计算gdop
end
