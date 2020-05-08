function energy_res_back=energy_distributionR(r,del,lamda,energy_res_back,z,depth)
%计算电子的能量沉积，传入参数为当前电子的所有r dE/dS, lamda,
%已有的能量沉积矩阵energy_res，z为电子纵坐标，depth为统计的深度最值

%统计出现次数，r对应的为横坐标，纵坐标加上值即可
    for i=1:length(r)
        temp=floor(r(i))+1;
        if temp>=1 && z(i)<depth && temp<2000
            energy_res_back(floor(r(i))+1)=energy_res_back(floor(r(i))+1)+del(i)*lamda(i);
        end
    end
end