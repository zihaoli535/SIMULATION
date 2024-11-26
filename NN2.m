function dx=NN2(t,x,xd)
c=50;e=0.05;a=0.5;
zeta=0.85;wn=80;
g11=1;g12=1/5;g21=1;g22=1/6.25;
L1=0.05;L2=0.05;ze1=2;ze2=1;b1=1;b2=1;o1=0.005;o2=0.005;
k11=5;k12=25;k21=6;k22=25;
m1=2;m2=2.5;j1=5;j2=6.25;b=0.5;
k=100;r=0.5;l=0.5;g=9.81;
mu=2;Tf=2;
y1r=0.5*sin(t)+0.5*sin(0.5*t);
dy1r=0.5*cos(t)+0.25*cos(0.5*t);
y2r=sin(0.5*t)*sin(t);
dy2r=0.5*cos(0.5*t)*sin(t)+sin(0.5*t)*cos(t);

cs1=-2;
cs2=-1;
cs3=0;
cs4=0.5;
cs5=1;
cs6=1.5;
X1=[x(1);x(2)];
X2=[x(3);x(4)];

phi11=exp(-norm(X1-cs1)^2/2)*ones(15,1);
phi21=exp(-norm(X2-cs1)^2/2)*ones(15,1);
phi12=exp(-norm(X1-cs2)^2/2)*ones(15,1);
phi22=exp(-norm(X2-cs2)^2/2)*ones(15,1);
phi13=exp(-norm(X1-cs3)^2/2)*ones(15,1);
phi23=exp(-norm(X2-cs3)^2/2)*ones(15,1);
phi14=exp(-norm(X1-cs4)^2/2)*ones(15,1);
phi24=exp(-norm(X2-cs4)^2/2)*ones(15,1);
phi15=exp(-norm(X1-cs5)^2/2)*ones(15,1);
phi25=exp(-norm(X2-cs5)^2/2)*ones(15,1);
phi16=exp(-norm(X1-cs6)^2/2)*ones(15,1);
phi26=exp(-norm(X2-cs6)^2/2)*ones(15,1);

if t<=Tf
     B=(1/t-1/Tf)^(2*mu)+e;
     dB=2*mu*(1/t-1/Tf)^(2*mu-1)*(-1/(t^2));
else
    B=e;
    dB=0;
end

z11=x(1)-y1r;
r1=c*(B^2-z11^2);
if r1>0&&r1<=a
    h1=1-(r1/a-1)^(2*mu);
    p1=1/h1-4*c*mu*(r1/a-1)^(2*mu-1)*z11^2/(a*h1^2);
    q1=4*c*mu*(r1/a-1)^(2*mu-1)*B*dB*z11/(a*h1^2);
else
    h1=1;
    p1=1;
    q1=0;
end

z21=x(3)-y2r;
r2=c*(B^2-z21^2);
if (r2>0)&&(r2<=a)
    h2=1-(r2/a-1)^(2*mu);
    p2=1/h2-4*c*mu*(r2/a-1)^(2*mu-1)*z21^2/(a*h2^2);
    q2=4*c*mu*(r2/a-1)^(2*mu-1)*B*dB*z21/(a*h2^2);
else
    h2=1;
    p2=1;
    q2=0;
end
f11=0; f12=((m1*g*r/j1-k*r^2/(4*j1))*sin(x(1))+k*r*(l-b)/(2*j1)); 
f21=0; f22=((m2*g*r/j2-k*r^2/(4*j2))*sin(x(3))+k*r*(l-b)/(2*j2));

z11=x(1)-y1r;  s11=z11/h1;  XI11=s11-x(9);
alp11=-(k11*s11+1.5*p1*XI11+p1*x(5)*f11-p1*dy1r+q1+p1^2*x(9))/(p1*g11);
dx13=wn*x(14);
dx14=-2*zeta*wn*x(14)-wn*(x(13)-alp11);
s12=x(2)-x(13);  XI12=s12-x(10);
u1=-(k12*s12+1.5*XI12+x(5)*f12+g11*s11*p1-dx13+(XI12*x(7)*norm(phi11))/(XI12^2+L1))/g12;

z21=x(3)-y2r;  s21=z21/h2;  XI21=s21-x(11);
alp21=-(k21*s21+1.5*p2*XI21+p2*x(6)*f21-p2*dy2r+q2+p2^2*x(11))/(p2*g21);
dx15=wn*x(16);
dx16=-2*zeta*wn*x(16)-wn*(x(15)-alp21);
s22=x(4)-x(15);  XI22=s22-x(12);
u2=-(k22*s22+1.5*XI22+x(6)*f22+g21*s21*p2-dx15+(XI22*x(8)*norm(phi21))/(XI22^2+L2))/g22;

dx1=x(2);%Example1
dx2=u1/j1+(m1*g*r/j1-k*r^2/(4*j1))*sin(x(1))+k*r*(l-b)/(2*j1)+k*r^2*sin(x(3))/(4*j1)+0.1*xd(1,1)*sin(xd(2,2))/j1+0.1*sin(t);
dx3=x(4);
dx4=u2/j2+(m2*g*r/j2-k*r^2/(4*j2))*sin(x(3))+k*r*(l-b)/(2*j2)+k*r^2*sin(x(1))/(4*j2)+0.1*xd(3,3)*cos(xd(4,4))/j1+0.1*cos(0.5*t);

dx5=p1*f11*XI11+f12*XI12-o1*x(5);
dx6=p2*f21*XI21+f22*XI22-o2*x(6);
dx7=ze1*XI12^2*norm(phi11)/(XI12^2+L1)-b1*ze1*x(7);
dx8=ze2*XI22^2*norm(phi21)/(XI22^2+L2)-b2*ze2*x(8);
dv11=-k11*x(9)+g11*p1*x(10)+g11*p1*(x(13)-alp11)-p1^2*x(9);
dv12=-k12*x(10)-g11*x(9);
dv21=-k21*x(11)+g21*p2*x(12)+g21*p2*(x(15)-alp21)-p2^2*x(11);
dv22=-k22*x(12)-g21*x(11);
dx72=(ze1*XI12^2*norm(phi12))/(XI12^2+L1)-b1*ze1*x(17);
dx82=(ze2*XI22^2*norm(phi22))/(XI22^2+L2)-b2*ze2*x(18);
dx73=(ze1*XI12^2*norm(phi13))/(XI12^2+L1)-b1*ze1*x(19);
dx83=(ze2*XI22^2*norm(phi23))/(XI22^2+L2)-b2*ze2*x(20);
dx74=(ze1*XI12^2*norm(phi14))/(XI12^2+L1)-b1*ze1*x(21);
dx84=(ze2*XI22^2*norm(phi24))/(XI22^2+L2)-b2*ze2*x(22);
dx75=(ze1*XI12^2*norm(phi15))/(XI12^2+L1)-b1*ze1*x(23);
dx85=(ze2*XI22^2*norm(phi25))/(XI22^2+L2)-b2*ze2*x(24);
dx76=(ze1*XI12^2*norm(phi16))/(XI12^2+L1)-b1*ze1*x(25);
dx86=(ze2*XI22^2*norm(phi26))/(XI22^2+L2)-b2*ze2*x(26);
dx=[dx1;dx2;dx3;dx4;dx5;dx6;dx7;dx8;dv11;dv12;dv21;dv22;dx13;dx14;dx15;dx16;dx72;dx82;dx73;dx83;dx74;dx84;dx75;dx85;dx76;dx86];
