clear
close all
tic

%定义一个cell，用来存储所有的参数
layer={};
%定义另一个cell，用来存储三维电子坐标以及能量信息
electron={};
hengzuobiao=150;
%用户选择输入方式
type=input('请选择键盘输入（0）或者文件读入（1）: ');
if type==0  
    %键盘输入
    %需要仿真的电子数目
    count_ele=input('请输入仿真的电子数目：');
    %绘制沿深度方向能量分布图像的范围（z的最大值）
    zmax=input('请输入绘制沿深度方向能量分布图像的范围（z的最大值）：');
    E=input('请输入电子能量初值');
    n=input('请输入层数: ');
    for i=1:n
        %Ci
        C_temp=input(['请输入第',num2str(i),'层材料的Ci: ']);   
        layer{i,1}=C_temp;
        %Zi
        Z_temp=input(['请输入第',num2str(i),'层材料的Zi: ']);   
        layer{i,2}=Z_temp;
        %Ai
        A_temp=input(['请输入第',num2str(i),'层材料的Ai: ']);
        layer{i,3}=A_temp;
        %rho
        rho_temp=input(['请输入第',num2str(i),'层材料的密度: ']);
        layer{i,4}=rho_temp;
        %depth
        zmax_temp=input(['请输入第',num2str(i),'层材料的厚度: ']);
        layer{i,5}=zmax_temp;
        %name
        name=input(['请输入第',num2str(i),'层材料的名称: '],'s');
        layer{i,6}=name;
    end
else
    %文件输入
    disp('文件内第一行为仿真电子数目，第二行为绘图时Z轴深度，第三行为电子能量E，之后每一行为参数，顺序为Ci Zi Ai rho depth');
    %以字符串方式读取键盘输入的文件名
    dataname=input('请输入文件名: ','s'); 
    file=fopen(dataname,'r');
    %需要仿真的电子数目
    count_ele=str2double(fgetl(file)); 
    %绘制沿深度方向能量分布图像的范围（z的最大值）
    zmax=str2double(fgetl(file)); 
    %获得电子能量E
    E=str2double(fgetl(file)); 
    %初始化n
    n=0; 
    %以是否到文件结尾为判断是否读取结束
    while ~feof(file)
        n=n+1;
        line_now=split(fgetl(file));    %将每一组数据切片
        layer{n,1}=eval(line_now{1});   %Ci
        layer{n,2}=eval(line_now{2});   %Zi
        layer{n,3}=eval(line_now{3});   %Ai
        layer{n,4}=eval(line_now{4});   %rho
        layer{n,5}=eval(line_now{5});   %depth       
        layer{n,6}=line_now{6};   %name
    end
    fclose(file);
end
%文件夹名称"depth_E_count_name1_name2"
filepath=[num2str(layer{1,5}),'_',num2str(E),'_',num2str(count_ele),'_',layer{1,6},'_',layer{2,6}];

%判断文件夹是否存在，如果不存在则新建此文件夹
if ~exist(['.\','result\',filepath],'dir')
    mkdir(['.\','result\',filepath]);
else
    str_temp=input('默认文件夹已存在，请输入文件夹名字末尾数字标识符，若不输入，则覆盖原文件：');
    filepath=[filepath,num2str(str_temp)];
    mkdir(['.\','result\',filepath]);
end

%将每一层的厚度存在一个名为depth的数组中
depth=zeros(1,n);
for n_layer=1:n
    depth(n_layer)=layer{n_layer,5};
end

%初始化能量沉积矩阵、背散射电子能量沉积矩阵、透射电子能量沉积矩阵
energy_res=zeros(1,floor(zmax)+1);
energy_res_back=zeros(1,hengzuobiao);
energy_res_forward=zeros(1,hengzuobiao);

%绘制水平线的横坐标范围
min_y=0;max_y=0;
min_x=0;max_x=0;
minx=[];maxx=[];miny=[];maxy=[];minz=[];maxz=[];
%背散射电子的数目
count_back=0;

%对每一个电子进行循环操作
for jjj=1:count_ele
    %每隔100个电子输出一次，让用户预估时间
    if mod(jjj,100)==0
        disp([num2str(jjj/100),' / ',num2str(count_ele/100)]);
    end
    
    %调用code_1计算并绘制电子的详细轨迹信息
    [energy_res_back,energy_res_forward,energy_res,miny_now,maxy_now,electron_temp,coor_temp,minx_now,maxx_now,count_back]=code_1(energy_res,zmax,E,layer,depth,energy_res_back,energy_res_forward,count_back,hengzuobiao);
    %更新最值
    if miny_now<min_y
        min_y=miny_now;
    end
    if maxy_now>max_y
        max_y=maxy_now;
    end
    
    if minx_now<min_x
        min_x=minx_now;
    end
    if maxx_now>max_x
        max_x=maxx_now;
    end
    
    %画能量分布二维图像所需的最值，以及每个电子的所有数据
    minx=[minx,coor_temp(1)];maxx=[maxx,coor_temp(2)];
    miny=[miny,coor_temp(3)];maxy=[maxy,coor_temp(4)];
    minz=[minz,coor_temp(5)];maxz=[maxz,coor_temp(6)];
    electron{jjj}=electron_temp;
end
%二维图像的最值
Min_x=min(minx);Max_x=max(maxx);
Min_y=min(miny);Max_y=max(maxy);
Min_z=min(minz);Max_z=max(maxz);

%最值有负数以及小数，变为整数间距，并且初始化二维矩阵
deltax=ceil(Max_x-Min_x+1);
deltay=ceil(Max_y-Min_y+1);
deltaz=ceil(Max_z-Min_z+1);
energy_distribution3=zeros(deltay,deltaz,1);

%对每一个电子进行操作
for k=1:count_ele
    %提取该电子数据  
    electron_temp=electron{k};
    for kk=1:length(electron_temp{1})
        %将二维矩阵对应的位置的数据与该电子能量沉积相加
        temp=energy_distribution3(ceil(electron_temp{2}(kk)+1-Min_y),ceil(electron_temp{3}(kk)+1-Min_z),1);
        temp=temp+electron_temp{4};
        energy_distribution3(ceil(electron_temp{2}(kk)+1-Min_y),ceil(electron_temp{3}(kk)+1-Min_z),1)=temp;   
    end

end
%变为正值方便绘制图像
energy_distribution3=-energy_distribution3;

%将数值变成最大值256，以绘制灰度图
energy_distribution3=256/max(max(energy_distribution3))*energy_distribution3;

%将每个小单元的能量累加，delta_y和delta_z为小单元的边长
delta_y=40;
delta_z=40;
%小单元格的数目
count_y=floor(deltay/delta_y); 
count_z=floor(deltaz/delta_z);
%对每一个小单元格进行操作，将内部元素的和作为内部所有元素新的值
for kx=1:count_y
    for kz=1:count_z
        energy_distribution3((kx-1)*delta_y+1:(kx)*delta_y,(kz-1)*delta_z+1:kz*delta_z)=ones(delta_y,delta_z).*sum(sum(energy_distribution3((kx-1)*delta_y+1:(kx)*delta_y,(kz-1)*delta_z+1:kz*delta_z))); 
    end
end

%除去边缘多出来的少量元素值
energy_distribution3=energy_distribution3(1:delta_y*count_y,1:delta_z*count_z);
%将最大值变为256，以绘制灰度图
energy_distribution3=256/max(max(energy_distribution3))*energy_distribution3;

%反转，使得背景为白色
energy_distribution3=-1.*energy_distribution3+256;
%energy_distribution3_back=-1.*energy_distribution3_back+256;
%最初存储的时候，横纵坐标反了，此处矩阵转置
energy_distribution3=energy_distribution3';

%绘图，电子能量沉积的二维灰度图
figure(2);imagesc(energy_distribution3);colormap(gray(256));title('energy distribution');

%指定向下绘制，故将z取反
energy_res=-energy_res;
energy_res_back=-energy_res_back;
energy_res_forward=-energy_res_forward;

%在电子轨迹的图像上绘制材料的分界线
for k=1:n-1
    if k==1
    %如果是第一条分界线，在此处分界线
        figure(5);line([min_x-200,max_x+200],[-(layer{k,5}-1),-(layer{k,5}+1)]);
        figure(1);[X,Y,z] = meshgrid(min_x-200:30:max_x+200,min_y-200:30:max_y+200,-(layer{k,5}));
    else
        figure(5);line([min_x-200,max_x+200],[-(layer{k,5}+layer{k-1,5}-1),-(layer{k,5}+layer{k-1,5}+1)]);
        figure(1);[X,Y,z] = meshgrid(min_x-200:30:max_x+200,min_y-200:30:max_y+200,-(layer{k,5}+layer{k-1,5}));
    end
end
mesh(X,Y,z)

%扩大沿Z轴以及R轴绘制的图像的横坐标间隔
energy_his=zeros(1,length(energy_res));
energy_his_back=zeros(1,length(energy_res_back));
energy_his_forward=zeros(1,length(energy_res_forward));

%指定沿深度方向能量沉积的横坐标间隔为50nm，沿R方向能量沉积的横坐标间距为20nm
delta_his=1;
delta_his_R=1;
%将每一个delta_his的能量沉积累加，作为这部分的值
for t=1:floor(zmax/delta_his)
    energy_his((t-1)*delta_his+1)=sum(energy_res((t-1)*delta_his+1:t*delta_his))/(count_ele*E);
end
for t=1:hengzuobiao/delta_his_R
    energy_his_back((t-1)*delta_his_R+1)=sum(energy_res_back((t-1)*delta_his_R+1:t*delta_his_R));
    energy_his_forward((t-1)*delta_his_R+1)=sum(energy_res_forward((t-1)*delta_his_R+1:t*delta_his_R));
end

for i=1:length(energy_his)
    figure(6);
    tempzong=linspace(0,energy_his(i),2);
    plot(i.*ones(length(tempzong)),tempzong,'color','b');hold on;
end
%画沿深度方向的能量沉积图
figure(3);plot([1:zmax+1],energy_his);title('energy distribution histogram');axis([-50,zmax+100,0,1.2*max(energy_his)]);xlabel('z/nm');
%画背散射电子、透射电子沿R方向的能量沉积图
figure(4);subplot(131);plot([1:hengzuobiao],energy_his_back);title('back');axis([-10,hengzuobiao,0,1.2*max(energy_his_back)]);xlabel('R/nm');
figure(4);subplot(132);plot([1:hengzuobiao],energy_his_forward);title('forward');axis([-10,hengzuobiao,0,1.2*max(energy_his_forward)]);xlabel('R/nm');
figure(4);subplot(133);plot([1:hengzuobiao],energy_his_back+energy_his_forward);title('all');axis([-10,hengzuobiao,0,1.2*max(energy_his_back+energy_his_forward)]);xlabel('R/nm');

%新建一个记录背散射系数的文本文件，将背散射系数记录在其中
backfile=fopen(['.\','result\',filepath,'\backscattered.txt'],'w');
fwrite(backfile,['Backscattering cofficient: ',num2str(count_back/count_ele)]);
fclose(backfile);
disp(['Backscattering cofficient: ',num2str(count_back/count_ele)]);

%将所有图像都保存到最初建立的文件夹之中，格式为png格式和fig格式两种
saveas(1,['.\','result\',filepath,'\trajectory3D.png']);
saveas(1,['.\','result\',filepath,'\trajectory3D.fig']);
saveas(2,['.\','result\',filepath,'\energy2D.png']);
saveas(2,['.\','result\',filepath,'\energy2D.fig']);
saveas(3,['.\','result\',filepath,'\energy1D.png']);
saveas(3,['.\','result\',filepath,'\energy1D.fig']);
saveas(5,['.\','result\',filepath,'\trajectory2D.png']);
saveas(5,['.\','result\',filepath,'\trajectory2D.fig']);
saveas(4,['.\','result\',filepath,'\energy_radius.png']);
saveas(4,['.\','result\',filepath,'\energy_radius.fig']);
saveas(6,['.\','result\',filepath,'\energy_z.png']);
saveas(6,['.\','result\',filepath,'\energy_z.fig']);
toc