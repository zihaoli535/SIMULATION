
function dx=Example2(t,x,xd)
c=100;e=0.05;a=0.5;
zeta=0.7;wn=270;
g11=1;g12=1;g21=1;g22=1;
L1=0.05;L2=0.05;ze1=10;ze2=10;b1=1;b2=1;
o1=0.01;o2=0.01;
k11=4;k12=18;k21=5;k22=20;
mu=2;Tf=1;
y1r=0.3*sin(t);
dy1r=0.3*cos(t);
y2r=0.3*cos(2*t);
dy2r=-0.6*sin(2*t);
cs=[-2 -1 1 2]';
X1=[x(1);x(2);x(5);x(7)];
X2=[x(3);x(4);x(6);x(8)];
phi1=exp(-norm(X1-cs)^2/2)*ones(15,1);
phi2=exp(-norm(X2-cs)^2/2)*ones(15,1);

if t<=Tf
     B=(1/t-1/Tf)^(2*mu)+e;
     dB=2*mu*(1/t-1/Tf)^(2*mu-1)*(-1/(t^2));
else
    B=e;
    dB=0;
end

z11=x(1)-y1r;
r1=c*(B^2-z11^2);
if (r1>0)&&(r1<=a)
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
f11=0; f12=(1-x(1)^2)*x(2)-x(1); f21=0; f22=(1-x(3)*x(4)^2);

z11=x(1)-y1r;  s11=z11/h1;  XI11=s11-x(9);
alp11=-(k11*s11+1.5*p1*XI11+p1*x(5)*f11-p1*dy1r+q1+p1^2*x(9))/(p1*g11);
dx13=wn*x(14);
dx14=-2*zeta*wn*x(14)-wn*(x(13)-alp11);
s12=x(2)-x(13);  XI12=s12-x(10);
u1=-(k12*s12+1.5*XI12+x(5)*f12+g11*s11*p1-dx13+(XI12*x(7)*norm(phi1))/(XI12^2+L1))/g12;

z21=x(3)-y2r;  s21=z21/h2;  XI21=s21-x(11);
alp21=-(k21*s21+1.5*p2*XI21+p2*x(6)*f21-p2*dy2r+q2+p2^2*x(11))/(p2*g21);
dx15=wn*x(16);
dx16=-2*zeta*wn*x(16)-wn*(x(15)-alp21);
s22=x(4)-x(15);  XI22=s22-x(12);
u2=-(k22*s22+1.5*XI22+x(6)*f22+g21*s21*p2-dx15+(XI22*x(8)*norm(phi2))/(XI22^2+L2))/g22;

dx1=x(2)+0.2*sin(x(3));
dx2=u1+(1-x(1)^2)*x(2)-x(1)+x(4)*sin(x(4))+xd(1,1)*sin(xd(2,2))+0.5*sin(t);
dx3=x(4)+0.2*sin(x(1));
dx4=u2+(1-x(3)*x(4)^2)+x(2)*cos(x(2))+xd(3,1)*cos(xd(4,2));

dx5=p1*f11*XI11+f12*XI12-o1*x(5);
dx6=p2*f21*XI21+f22*XI22-o2*x(6);
dx7=(ze1*XI12^2*norm(phi1))/(XI12^2+L1)-b1*ze1*x(7);
dx8=(ze2*XI22^2*norm(phi2))/(XI22^2+L2)-b2*ze2*x(8);
dv11=-k11*x(9)+g11*p1*x(10)+g11*p1*(x(13)-alp11)-p1^2*x(9);
dv12=-k12*x(10)-g11*x(9);
dv21=-k21*x(11)+g21*p2*x(12)+g21*p2*(x(15)-alp21)-p2^2*x(11);
dv22=-k22*x(12)-g21*x(11);

dx=[dx1;dx2;dx3;dx4;dx5;dx6;dx7;dx8;dv11;dv12;dv21;dv22;dx13;dx14;dx15;dx16];
