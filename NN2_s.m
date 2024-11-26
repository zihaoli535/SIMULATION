clc;close all;
sol=dde23(@NN2,[0.2 0.3 0.2 0.3],[1;1.3;0.5;1.3;0;0;0;0;0;0;0;0;-6.2;0;-3.45;0;0;0;0;0;0;0;0;0;0;0],[0 10]);

t=sol.x;
x=sol.y;
c=50;e=0.05;a=0.5;
zeta=0.85;wn=80;
g11=1;g12=1/5;g21=1;g22=1/6.25;
L1=0.05;L2=0.05;
k11=5;k12=25;k21=6;k22=25;
m1=2;m2=2.5;j1=5;j2=6.25;b=0.5;
k=100;r=0.5;l=0.5;g=9.81;
mu=2;Tf=2;
y1r=0.5*sin(t)+0.5*sin(0.5*t);
dy1r=0.5*cos(t)+0.25*cos(0.5*t);
y2r=sin(0.5*t).*sin(t);
dy2r=0.5*cos(0.5*t).*sin(t)+sin(0.5*t).*cos(t);

Bv = ((1./(t+0.001) - 1./2.001).^4 + 0.05) .* (t <= 2) + 0.05 .* (t > 2);
cs=[-2 -1 1 2]';
X1=[x(1,:);x(2,:);x(5,:);x(7,:)];
X2=[x(3,:);x(4,:);x(6,:);x(8,:)];
phi1=exp(-norm(X1-cs).^2./2).*ones(15,1);
phi2=exp(-norm(X2-cs).^2./2).*ones(15,1);

u1 = zeros(length(t),1);
u2 = zeros(length(t),1);
for i = 1:length(t)
if t<=Tf
     B=(1./t-1./Tf).^(2.*mu)+e;
     dB=2.*mu.*(1./t-1./Tf).^(2.*mu-1).*(-1./(t.^2));
else
    B=e;
    dB=0;
end
 
z11=x(1,:)-y1r;
r1=c.*(B.^2-z11.^2);
if  (r1>0)&(r1<=a)
    h1=1-(r1./a-1).^(2.*mu);
    p1=1./h1-4.*c.*mu.*(r1./a-1).^(2.*mu-1).*z11.^2./(a.*h1.^2);
    q1=4.*c.*mu.*(r1./a-1).^(2.*mu-1).*B.*dB.*z11./(a.*h1.^2);
else
    h1=1;
    p1=1;
    q1=0;
end
 
z21=x(3,:)-y2r;
r2=c.*(B.^2-z21.^2);
if  (r2>0)&(r2<=a)
    h2=1-(r2./a-1).^(2.*mu);
    p2=1./h2-4.*c.*mu.*(r2./a-1).^(2.*mu-1).*z21.^2./(a.*h2.^2);
    q2=4.*c.*mu.*(r2./a-1).^(2.*mu-1).*B.*dB.*z21./(a.*h2.^2);
else
    h2=1;
    p2=1;
    q2=0;
end
f11=0; f12=((m1*g*r/j1-k*r^2/(4*j1))*sin(x(1,:))+k*r*(l-b)/(2*j1));
f21=0; f22=((m2*g*r/j2-k*r^2/(4*j2))*sin(x(3,:))+k*r*(l-b)/(2*j2));

z11=x(1,:)-y1r;  s11=z11./h1;  XI11=s11-x(9,:);
alp11=-(k11.*s11+1.5.*p1.*XI11+p1.*x(5,:).*f11-p1.*dy1r+q1+p1^2*x(9,:))./(p1.*g11);

dx13=wn*x(14,:);
dx14=-2*zeta*wn*x(14,:)-wn*(x(13,:)-alp11);
s12=x(2,:)-x(13,:);  XI12=s12-x(10,:);
u1=-(k12.*s12+1.5.*XI12+x(5,:).*f12+g11.*s11.*p1-dx13+(XI12.*x(7,:).*norm(phi1))./(XI12.^2+L1))./g12;

z21=x(3,:)-y2r;  s21=z21./h2;  XI21=s21-x(11,:);
alp21=-(k21.*s21+1.5.*p2.*XI21+p2.*x(6,:).*f21-p2.*dy2r+q2+p2^2*x(11,:))./(p2.*g21);

dx15=wn*x(16,:);
dx16=-2*zeta*wn*x(16,:)-wn*(x(15,:)-alp21);
s22=x(4,:)-x(15,:);  XI22=s22-x(12,:);
u2=-(k22.*s22+1.5.*XI22+x(6,:).*f22+g21.*s21.*p2-dx15+XI22.*x(8,:).*norm(phi2)./(XI22.^2+L2))./g22;
end
figure(1)
plot(t,x(7,:));
hold on
plot(t,x(17,:));
hold on
plot(t,x(19,:));
hold on
plot(t,x(21,:));
hold on
plot(t,x(23,:));
hold on
plot(t,x(25,:));
legend('$\hat\eta_{1,1}$','$\hat\eta_{1,2}$','$\hat\eta_{1,3}$','$\hat\eta_{1,4}$','$\hat\eta_{1,5}$','$\hat\eta_{1,6}$','interpreter','latex','fontsize',18)
ylim([0 0.15]);
xlabel('Time(sec)')

figure(2)
plot(t,x(8,:));
hold on
plot(t,x(18,:));
hold on
plot(t,x(20,:));
hold on
plot(t,x(22,:));
hold on
plot(t,x(24,:));
hold on
plot(t,x(26,:));
legend('$\hat\eta_{2,1}$','$\hat\eta_{2,2}$','$\hat\eta_{2,3}$','$\hat\eta_{2,4}$','$\hat\eta_{2,5}$','$\hat\eta_{2,6}$','interpreter','latex','fontsize',18)
ylim([0 0.15]);
xlabel('Time(sec)')

