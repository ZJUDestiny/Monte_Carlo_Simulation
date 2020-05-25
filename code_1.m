function [energy_res_back,energy_res_forward,energy_res,miny,maxy,electron_temp,coor_temp,minx,maxx,count_back]=code_1(energy_res,zmax,E,layer,depth_all,energy_res_back,energy_res_forward,count_back,hengzuobiao)
%输入参数：已有的能量分布矩阵energy_res，z的最大值zmax，电子能量初始能量E，包含所有参数的layer，所有层的深度depth_all
%背散射电子的能量沉积矩阵energy_res_back，透射电子的能量沉积矩阵energy_res_forward，背散射电子的数目count_back

%输出参数：加上此电子的能量分布矩阵energy_res_back, energy_res_forward, energy_res，
%此电子的坐标的最小值minx miny、最大值maxx maxy，当前电子的所有数据electron_temp，坐标数据electron_temp，
%背散射电子数目count_back

    %定义是否背散射
    is_back=0;
    
    %上一次的depth，用于累加depth，计算是否越过边界
    depth_last=0;
    %初始化能量分布的所需矩阵
    %每一个电子的x,y,z矩阵
    x_x=[];
    y_x=[];
    z_x=[];
    %每一个电子对应的z的dE/dS
    del_energy=[];
    %每一个电子对应的z的lamda
    lamda_energy=[];
    
    %初始化电子最终位置的矩阵
    finalx=[];
    finaly=[];
    finalz=[];

    %阿伏伽德罗常数
    NA=6.02*10^23; 

    %计算初始坐标
    %电子束半径1nm，并把单位换算成纳米
    Re=1; 
    Re=Re*10^-7; 
    %两个零到一的随机数
    w1=rand();w2=rand(); 
    x0=Re*sqrt(-log(w1))*cos(2*pi*w2);
    y0=Re*sqrt(-log(w1))*sin(2*pi*w2);
    z0=0;
    
    %将当前坐标存储到电子坐标的矩阵之中，乘上10^7变为纳米单位
    x_now=x0;y_now=y0;z_now=z0;
    finalx=[finalx,x_now*10^7];finaly=[finaly,y_now*10^7];finalz=[finalz,z_now*10^7];

    %获取第一种材料的参数值
    count=1;
    Ci=layer{count,1};
    Zi=layer{count,2};
    Ai=layer{count,3};
    rho=layer{count,4};
    depth=depth_last+layer{1,5};
    
    %depth_add，获得当前层深度的绝对值
    depth_add=zeros(1,length(depth_all)+1);
    for n_depth=1:length(depth_add)-1
        depth_add(n_depth+1)=sum(depth_all(1:n_depth));
    end
    
    %判断是否碰撞结束，条件为E>0.2keV
    while E>0.2
        %判断是否进入了下一种材料之中,用电子坐标是否超过已有的所有depth之和，而不是采用当前材料的depth
        for layer_count=1:length(depth_add)-1
            %判断电子进入了哪一层的材料之中，并更新新的一层材料参数
            if z_now*10^7>=depth_add(layer_count) && z_now*10^7<depth_add(layer_count+1)
                count=layer_count;
                Ci=layer{count,1};
                Zi=layer{count,2};
                Ai=layer{count,3};
                rho=layer{count,4};
                depth_last=depth;
                depth=depth_last+layer{count,5};
                break;
            end
        end   
        
        %屏蔽参数，E-keV
        beta_i=5.43*10^-3.*Zi.^(2/3)./E; 
        %得到屏蔽的Rutherford弹性散射截面，单位cm2
        sigma_i=5.21*10^-21.*Zi.*(Zi+1)./E^2.*((E+511)./(E+1022))^2*4*pi./beta_i./(beta_i+1); 
        %多原子散射截面
        ni_sigmai=Ci.*rho.*NA./Ai.*sigma_i;  
        %平均自由程
        lamda_ba=1/sum(ni_sigmai);
        %与每种原子碰撞的概率
        P=ni_sigmai./sum(ni_sigmai);  
        %判断在当前物质中和谁碰撞
        %用来判断与哪个原子碰撞的随机数R0
        R0=rand();
        %获取当前物质原子种类数目
        n_atom=length(Ci);
        %如果小于第一种的碰撞概率，i=1
        if 0<=R0 && R0<=P(1)
            i=1;
        else
        %否则进行判断，当大于已有的原子概率和，小于加上下一种原子的概率和，和下一种原子碰撞，i=nn+1
            for nn=1:n_atom-1
                if R0>=sum(P(1:nn)) && R0<sum(P(1:nn+1))
                    i=nn+1;
                end
            end
        end
        %获取碰撞的原子的屏蔽系数，并且计算自由程
        beta=beta_i(i);     
        lamda=-lamda_ba*log(rand());
        %散射角theta、方位角fai
        [theta,fai]=cal_angle(beta);
        
        %计算能量损耗
        %调用函数计算Ji
        Ji=cal_J(Zi);
        
        %计算化合物的平均原子序数、平均原子量、平均电离能
        %将化合物等效成Z_ba A_ba J_ba的单质
        Z_ba=sum(Ci.*Zi./Ai)/sum(Ci./Ai);
        A_ba=1/sum(Ci./Ai);
        J_ba=exp(sum(Ci.*Zi./Ai.*log(Ji))/sum(Ci.*Zi./Ai));
      
        %计算dE/dS，以及后续能量，其中包括了修正系数
        del=-7.85*10^4*rho*Z_ba/A_ba/E*log(1.166*E/J_ba+0.0063*Z_ba+0.7334);
        
        %计算坐标转换
        %判断如果转换矩阵Q还不存在，那么初始化为单位矩阵
        if ~exist('Q')
            Q=eye(3);
        else
        %如果存在，那么调用函数计算Qn+1,n，并与传入的Qn,0做矩阵乘法，求出Qn+1,0
            Q=cal_Q(Q,theta,fai);
        end
        
        %将x,y,z坐标、dE/dS以及lamda均保存下来以计算能量分布
        x_x=[x_x,x_now*10^7];
        y_x=[y_x,y_now*10^7];
        z_x=[z_x,z_now*10^7];
        
        %保存上次的电子坐标
        x_last=x_now;y_last=y_now;z_last=z_now;
        %计算碰撞后的电子位置
        [x_now,y_now,z_now]=cal_newcoor(Q,lamda,theta,fai,x_now,y_now,z_now); 
        
        
        %首先判断层数，如果跨层则需要 计算全新参数
        for layer_count_2=1:length(depth_add)-1
            if z_now*10^7>=depth_add(layer_count_2) && z_now*10^7<depth_add(layer_count_2+1)
                %判断层数
                count_2=layer_count_2;
                if count_2~=count
                    %出现跨层现象
                    %新的一层的数据
                    Ci_2=layer{count_2,1};
                    Zi_2=layer{count_2,2};
                    Ai_2=layer{count_2,3};
                    rho_2=layer{count_2,4};
                    %屏蔽参数，E-keV
                    beta_i_2=5.43*10^-3.*Zi_2.^(2/3)./E; 
                    %得到屏蔽的Rutherford弹性散射截面，单位cm2
                    sigma_i_2=5.21*10^-21.*Zi_2.*(Zi_2+1)./E^2.*((E+511)./(E+1022))^2*4*pi./beta_i_2./(beta_i_2+1); 
                    %多原子散射截面
                    ni_sigmai_2=Ci_2.*rho_2.*NA./Ai_2.*sigma_i_2;  
                    %平均自由程
                    lamda_ba_2=1/sum(ni_sigmai_2);
                    %散射步长
                    
                    Ji_2=cal_J(Zi_2);
                    %计算化合物的平均原子序数、平均原子量、平均电离能
                    %将化合物等效成Z_ba A_ba J_ba的单质
                    Z_ba_2=sum(Ci_2.*Zi_2./Ai_2)/sum(Ci_2./Ai_2);
                    A_ba_2=1/sum(Ci_2./Ai_2);
                    J_ba_2=exp(sum(Ci_2.*Zi_2./Ai_2.*log(Ji_2))/sum(Ci_2.*Zi_2./Ai_2));
                    %计算校正后的dE/dS，以及后续能量
                    del_2=-7.85*10^4*rho_2*Z_ba_2/A_ba_2/E*log(1.166*E/J_ba_2+0.0063*Z_ba_2+0.7334);
                    
                    %算出原来的两点之间距离，与lamda应该一致
                    l=sqrt(((x_last-x_now))^2+((y_last-y_now))^2+((z_last-z_now))^2);
                    

                    %算出接触面以前的长度，根据公式计算后一段长度
                    l1=(z_last-depth_add(max(count_2,count))*10^-7)/(z_last-z_now)*l;
                    lamda_2=(lamda-l1)*lamda_ba_2/lamda_ba;
                    %计算修正后的散射步长lamda
                    lamda= l1+lamda_2;
                    %计算修正后的能量损失率
                    del=(del*l1+del_2*lamda_2)/lamda;
                    
                    %计算修正后的电子散射终点的坐标
                    x_now=(x_now-x_last)*lamda/l+x_last;
                    y_now=(y_now-y_last)*lamda/l+y_last;
                    z_now=(z_now-z_last)*lamda/l+z_last;
                        
                end
                break;
            end
        end  
        
        %将电子能量计算出来，把del和lamda存起来便于画图
        E=E+del*lamda;
        del_energy=[del_energy,del];
        lamda_energy=[lamda_energy,lamda];
      
        %设置了边界条件
        if z_now*10^7>0
            %如果还在固体之中，则存下坐标
            finalx=[finalx,x_now*10^7];finaly=[finaly,y_now*10^7];finalz=[finalz,z_now*10^7];
        else
            %如果已经逸出固体表面，则把E设为零，并将其标记为背散射电子
            E=0;
            is_back=1;
        end
    end
    %将电子轨迹绘图
    figure(1);title('Electronic trajectory3D');xlabel('x/nm');ylabel('y/nm');zlabel('z/nm');
    %判断是否背散射电子，用不同颜色标注
    if is_back==1
        plot3(finalx,finaly,-finalz,'r');hold on
    else
        plot3(finalx,finaly,-finalz,'b');hold on
    end   
    figure(5);title('Electronic trajectory');xlabel('x/nm');ylabel('z/nm');
    if is_back==1
        count_back=count_back+1;
        plot(finaly,-finalz,'r');hold on
    else
        plot(finaly,-finalz,'b ');hold on
    end 
    %返回此电子的横坐标最值
    miny=min(finaly);maxy=max(finaly);
    minx=min(finalx);maxx=max(finalx);

    %将最后一次散射后的电子坐标以及能量损失率、能量损失皆保存下来
    x_x=[x_x,x_now*10^7];
    y_x=[y_x,y_now*10^7];
    z_x=[z_x,z_now*10^7];
    del_energy=[del_energy,E];
    lamda_energy=[lamda_energy,-1];
    %计算能量分布
    energy_res=energy_distribution(z_x,del_energy,lamda_energy,zmax,energy_res);
    %画横坐标是r的时候的能量图
    r_x=sqrt(x_x.^2+y_x.^2);
    if is_back==1
        %背散射电子
        energy_res_back=energy_distributionR(r_x,del_energy,lamda_energy,energy_res_back,z_x,layer{1,5},hengzuobiao);
    else
        %前散射电子
        energy_res_forward=energy_distributionR(r_x,del_energy,lamda_energy,energy_res_forward,z_x,layer{1,5},hengzuobiao);
    end
    %将此次电子散射的数据均保存到一个cell中，第一个数据为所有散射的x坐标，第二个为y，第三个为z，第四个为每次散射损失的能量
    electron_temp={};
    electron_temp{1}=finalx;
    electron_temp{2}=finaly;
    electron_temp{3}=finalz;
    electron_temp{4}=del.*lamda;
    %将电子坐标最值返回以方便绘图
    coor_temp=[min(finalx),max(finalx),min(finaly),max(finaly),min(finalz),max(finalz)];
end

