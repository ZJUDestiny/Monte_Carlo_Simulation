function energy_res=energy_distribution(z,del,lamda,zmax,energy_res)
%计算电子的能量沉积，传入参数为当前电子的所有z dE/dS, lamda, 沉积的最大zmax，
%已有的能量沉积矩阵energy_res，z为电子纵坐标，depth为统计的深度最值

%统计出现次数，z对应的为横坐标，纵坐标加上值即可
    for i=1:length(z)
        temp=floor(z(i))+1;
        if temp<=zmax && temp>1
            energy_res(floor(z(i))+1)=energy_res(floor(z(i))+1)+del(i)*lamda(i);
        end
    end
end