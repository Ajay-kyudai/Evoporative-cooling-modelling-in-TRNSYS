clc
clear all
% constants for beta
A1=2.905392768;
A2=-3.204944487;
A3=1.320700525;
phi_max=1.1;
% Theoritical simulation
n=40; 
T=341940; 
dt=5;
m=T/dt;
m1=m/12;
m2=5700;
m3=1141;
t1=linspace(0,T,m);
a=0.2; % m, Thickness of roof
dx=zeros(n,1);
Tini=296;          % K        
L=2257000;         % Latent heat of vaporization of water, J/Kg 
alpha = 0.6;      % absorptivity of concrete
alpharoom=9.67;       % inside heat transfer coefficient W/m2/K
sigma=5.67e-08;    % stefan boltzman constant
emissivity=0.94;   % emissivity of concrete
lewis=0.83;           % lewis number
rhoda=1.2041;      % density of dry air
t=0; % initial time
phiini=0.06; % initial water content
phi_tk=phiini;
Troom=24; % degree C
% data reading
sw=zeros(m2,1);
Tout=zeros(m2,1);
AH=zeros(m2,1);
Pr=zeros(m2,1);
v=zeros(m2,1);
P=zeros(m2,1);
w=xlsread('Jodhpur Jun minutewise data.xlsx');
for i=1:m2
    sw(i)=w(i,2);
end
for i=1:m2
    Tout(i)=w(i,3);
end
for i=1:m2
    AH(i)=w(i,4);
end
for i=1:m2
    Pr(i)=w(i,7);
end
for i=1:m2
    v(i)=w(i,5);
end
for i=1:m2
    P(i)=w(i,6);
end
% Interpolated data
sw01=zeros(m,1);
Tout01=zeros(m,1);
AH01=zeros(m,1);
Pr01=zeros(m,1);
v01=zeros(m,1);
P01=zeros(m,1);
BETA=zeros(m,1);
phi= zeros(m,1);
R=zeros(42,m);
EVflux=zeros(m,1);
netradiation=zeros(m,1);
CV=zeros(m,1);
CD=zeros(m,1);
Balance=zeros(m,1);
Thermalload=zeros(m,1);
watermass=zeros(m,1);
Xsat=zeros(m,1);
KX=zeros(m,1);
for i=1:m2-1
    k=12*(i-1)+1;
    for j=k:12*i
        r1=(sw(i+1,:)-sw(i,:))/11;
        sw01(j,:)=sw(i,:)+(j-12*(i-1)-1)*r1;
    end
end
for i=1:m2-1
    k=12*(i-1)+1;
    for j=k:12*i
        r2=(Tout(i+1,:)-Tout(i,:))/11;
        Tout01(j,:)=Tout(i,:)+(j-12*(i-1)-1)*r2;
    end
end
for i=1:m2-1
    k=12*(i-1)+1;
    for j=k:12*i
        r3=(AH(i+1,:)-AH(i,:))/11;
        AH01(j,:)=AH(i,:)+(j-12*(i-1)-1)*r3;
    end
end
for i=1:m2-1
    k=12*(i-1)+1;
    for j=k:12*i
        r4=(Pr(i+1,:)-Pr(i,:))/11;
        Pr01(j,:)=Pr(i,:)+(j-12*(i-1)-1)*r4;
    end
end
for i=1:m2-1
    k=12*(i-1)+1;
    for j=k:12*i
        r5=(v(i+1,:)-v(i,:))/11;
        v01(j,:)=v(i,:)+(j-12*(i-1)-1)*r5;
    end
end
for i=1:m2-1
    k=12*(i-1)+1;
    for j=k:12*i
        r6=(P(i+1,:)-P(i,:))/11;
        P01(j,:)=P(i,:)+(j-12*(i-1)-1)*r6;
    end
end
% Initialization for 1 minute interval
R1=zeros(42,m1);
EVflux1=zeros(m1,1);
netradiation1=zeros(m1,1);
CV1=zeros(m1,1);
CD1=zeros(m1,1);
Balance1=zeros(m1,1);
Thermalload1=zeros(m1,1);
watermass1=zeros(m1,1);
phi1=zeros(m1,1);
sw1=zeros(m1,1);
AH1=zeros(m1,1);
P1=zeros(m1,1);
V1=zeros(m1,1);
Tout1=zeros(m1,1);
Pr1=zeros(m1,1);
Xsat1=zeros(m1,1);
BETA1=zeros(m1,1);
KX1=zeros(m1,1);
% calculation of averaged values at each 5 minute interval
netradiation2=zeros(m3,1);
CV2=zeros(m3,1);
CD2=zeros(m3,1);
EVflux2=zeros(m3,1);
Ttop1=zeros(m3,1);
Tbottom1=zeros(m3,1);
Thermalload2=zeros(m3,1);
Pr2=zeros(m3,1);
BETA2=zeros(m3,1);
phi2=zeros(m3,1);
 % creation of initial temperature matrix U0
for i=1:n
    dx(i)=0.2/n;
end 
U_tk=zeros(n+2,1);
    for i=1:n+2
        U_tk(i)=Tini;
    end
U=zeros(n+2,m);
U(:,1)=U_tk;
% loop starts with time step of 5s
 for l=1:m
     lambda=zeros(n,1);
    rho=zeros(n+2,1);
    cp=zeros(n+2,1);
    
   for i=1:n
        lambda(i)= 1.58;
        rho(i)=2248;
        cp(i)=880;
   end
     
        alphao=3.96*v01(l)+6.42;
    
    % creation of heat capacity matrix M
    M=zeros(n+2,n+2);
    for i=2:n+1
        for j=2:n+1
            if i==j
                M(i,j)=rho(i-1)*cp(i-1)*dx(i-1);
            end
        end
    end
    
    %creation of C
    c=zeros(n+1,1);
    c(1)=(lambda(1)/(dx(1)/2));
    for i=2:n
        c(i)= (1/((dx(i-1)/(2*lambda(i-1)))+(dx(i)/(2*lambda(i)))));
    end
    c(n+1)=(2*lambda(n)/dx(n));
    
    % creation of heat transfer coefficient matrix S 
   S= zeros(n+2,n+2);
   for i=2:n+1
       for j=2:n+1
           if i==j
               S(i,j)=-c(i-1)-c(i);
           end
       end
   end
   
   S(1,1)= -c(1)-alphao;
   S(n+2,n+2)= -c(n+1)-alpharoom;
   for i=1:n+1
       for j=1:n+2
           if j==i+1
               S(i,j)=c(i);
           end
       end
   end
   for i=2:n+2
       for j=1:n+1
           if j==i-1
               S(i,j)=c(i-1);
           end
       end
   end
   M1=(1/dt)*M;
   A=M1-S;
     
    % creation of convection boundary condition matrix E
    D= zeros(n+2,2);
    D(1,1)= alphao;
    D(n+2,2)= alpharoom;
    Theta=zeros(2,1);
    Theta(1)=Tout01(l)+273.15;
    Theta(2)=Troom+273.15;
    E=D*Theta;
   tvec=zeros(m,1);
   
   tvec(1)=t;
   F=zeros(n+2,1);
   
   % complete definition of grid
   
   x=zeros(n+2,1);
   x(1)=0.0;
   x(2)=dx(1)/2;
   for i=2:40
      x(i+1)=x(i)+0.5*(dx(i-1)+dx(i));
    end
 
   x(n+2)=a;
   % creation of main loop
        t=t+dt;
        Pwss=6.116441*10^(7.591386*(U_tk(1)-273.15)/((U_tk(1)-273.15)+240.7263));
        Xss=0.622*Pwss/(P01(l)-Pwss);
        Pw=P01(l)*AH01(l)/(AH01(l)+0.622);
        f=Pw*10^2/(13.6*9.81);
        epsky= 0.526+0.076*sqrt(f);
        Isky=epsky*sigma*(Tout01(l)+273.15)^4;
        Rnet=alpha*sw01(l)+emissivity*Isky-emissivity*sigma*(U_tk(1))^4;
        Cpm=1000+AH01(l)*2100;
        rhoma=rhoda*(1+AH01(l))/(1+1.609*AH01(l));
        Kx=alphao/(lewis*rhoma*Cpm);
        phi_tk_1=phi_tk+dt*(Pr01(l)-((A1*(phi_tk/phi_max)^3+A2*(phi_tk/phi_max)^2+A3*(phi_tk/phi_max))*Kx*(Xss-AH01(l))));
        if phi_tk_1<= phi_max
            phi(l)=phi_tk_1;
            phi_tk=phi_tk_1;
        else
            phi_tk_1=phi_max;
            phi(l)=phi_tk_1;
            phi_tk=phi_tk_1;
        end 
        beta=A1*(phi_tk/phi_max)^3+A2*(phi_tk/phi_max)^2+A3*(phi_tk/phi_max);
        EV=beta*Kx*(Xss-AH01(l));
        F(1,1)=Rnet-L*EV; 
        Z=(M1*U_tk)+E+F;
        U_tk_1=A\Z;
        U(:,l)=U_tk_1;
        U_tk=U_tk_1;
        tvec(l)=t;
  KX(l,:)=Kx;
 BETA(l,:)=beta;
  phi(l,1)=phi_tk;
   U(:,l)=U_tk;
   R(:,l)=U(:,l);
   EVflux(l,:)=L*EV;
   Xsat(l,:)=Xss;
   netradiation(l,:)=Rnet;
   CV(l,:)=alphao*((Tout01(l)+273.15)-U_tk(1));
   CD(l,:)=c(1)*(U_tk(2)-U_tk(1));
   Balance(l,:)=CV(l,:)+CD(l,:)+netradiation(l,:)- EVflux(l,:);
   Thermalload(l,:)=alpharoom*(U_tk(n+2)-(Troom+273.15));
   if EV>0
       watermass(l,:)=EV;
   end
 end
 % to get output at 1 minute interval
 for k=1
   z=k;
    R1(:,k)=R(:,z);
    EVflux1(k,:)=EVflux(z,:);
    CV1(k,:)=CV(z,:);
    CD1(k,:)=CD(z,:);
    netradiation1(k,:)=netradiation(z,:);
    Balance1(k,:)=Balance(z,:);
    Thermalload1(k,:)=Thermalload(z,:);
    watermass1(k,:)=watermass(z,:);
    phi1(k,:)=phi(z,:);
    sw1(k,:)=sw01(z,:);
    AH1(k,:)=AH01(z,:);
    V1(k,:)=v01(z,:);
    Tout1(k,:)=Tout01(z,:);
    Pr1(k,:)=Pr01(z,:);
   P1(k,:)=P01(z,:);
   Xsat1(k,:)=Xsat(z,:);
   BETA1(k,:)=BETA(z,:);
   KX1(k,:)=KX(z,:);
 end
for k=2:m1+1
    z=(m/m1)*(k-1);
    R1(:,k)=R(:,z);
    EVflux1(k,:)=EVflux(z,:);
    CV1(k,:)=CV(z,:);
    CD1(k,:)=CD(z,:);
    netradiation1(k,:)=netradiation(z,:);
    Balance1(k,:)=Balance(z,:);
    Thermalload1(k,:)=Thermalload(z,:);
    watermass1(k,:)=watermass(z,:);
    phi1(k,:)=phi(z,:);
    sw1(k,:)=sw01(z,:);
    AH1(k,:)=AH01(z,:);
    V1(k,:)=v01(z,:);
    Tout1(k,:)=Tout01(z,:);
    Pr1(k,:)=Pr01(z,:);
   P1(k,:)=P01(z,:);
   Xsat1(k,:)=Xsat(z,:);
   BETA1(k,:)=BETA(z,:);
   KX1(k,:)=KX(z,:);
end 
R2=R1';
Ttop=R2(:,1);
Tbottom=R2(:,42);
% 5 minute avg result

    EVflux2(1)=EVflux1(1);
    CV2(1)=CV1(1);
    CD2(1)=CD1(1);
    netradiation2(1)=netradiation1(1);
    Thermalload2(1)=Thermalload1(1);
    Ttop1(1)=Ttop(1);
    Tbottom1(1)=Tbottom(1);
    Pr2(1)=Pr1(1);
    BETA2(1)=BETA1(1);
    phi2(1)=phi1(1);

 for k=2:m3
    z1=5*(k-2)+1;
    z2=5*(k-1);
    EVflux2(k,:)=mean(EVflux1(z1:z2));
    CV2(k,:)=mean(CV1(z1:z2));
    CD2(k,:)=mean(CD1(z1:z2));
    netradiation2(k,:)=mean(netradiation1(z1:z2));
    Thermalload2(k,:)=mean(Thermalload1(z1:z2));
    Ttop1(k,:)=mean(Ttop(z1:z2));
    Tbottom1(k,:)=mean(Tbottom(z1:z2));
    Pr2(k,:)=mean(Pr1(z1:z2));
    BETA2(k,:)=mean(BETA1(z1:z2));
    phi2(k,:)=mean(phi1(z1:z2));
end 

figure(2)
plot (t1,phi)
xlabel('time');
ylabel('phi')
figure(3)
plot(x,U(:,m))
xlabel('roof suface depth')
ylabel('Temperature(K)')
