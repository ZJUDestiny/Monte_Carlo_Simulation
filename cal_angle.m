function [theta,fai]=cal_angle(beta)
    %散射角theta，方位角fai的计算函数
    R=rand();
    costheta=1-2*beta*R/(1+beta-R);
    %防止cos值小于-1，程序错误
    if costheta<-1
        costheta=-1;
    end
    theta=acos(costheta); %弧度制
    R=rand();
    fai=2*pi*R;
end