function [x_now,y_now,z_now]=cal_newcoor(Q,lamda,theta,fai,x0,y0,z0)
%����lamda theta fai���x y z
    x_localcoor=lamda*sin(theta)*cos(fai);
    y_localcoor=lamda*sin(theta)*sin(fai);
    z_localcoor=lamda*cos(theta);
%���x y z������ 
    delta=Q*[x_localcoor;y_localcoor;z_localcoor];
    deltax0=delta(1);deltay0=delta(2);deltaz0=delta(3);
%����һ��������ӣ�����µĵ�������
    x_now=deltax0+x0;
    y_now=deltay0+y0;
    z_now=deltaz0+z0;
end