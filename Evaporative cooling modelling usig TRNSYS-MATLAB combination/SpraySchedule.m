
% trnTime (1x1)        : simulation time 
% trnInfo (15x1)       : TRNSYS info array
% trnInputs (nIx1)     : TRNSYS inputs 
% trnStartTime (1x1)   : TRNSYS Simulation Start time
% trnStopTime (1x1)    : TRNSYS Simulation Stop time
% trnTimeStep (1x1)    : TRNSYS Simulation time step
% mFileErrorCode (1x1) : Error code for this m-file. It is set to 1 by TRNSYS and the m-file should set it to 0 at the
%                        end to indicate that the call was successful. Any non-zero value will stop the simulation
% trnOutputs (nOx1)    : TRNSYS outputs 
%trnInputs             : TRNSYS Inputs

% TRNSYS sets mFileErrorCode = 1 at the beginning of the M-File for error detection
% This file increments mFileErrorCode at different places. If an error occurs in the m-file the last succesful step will
% be indicated by mFileErrorCode, which is displayed in the TRNSYS error message
% At the very end, the m-file sets mFileErrorCode to 0 to indicate that everything was OK

mFileErrorCode = 100    % Beginning of the m-file 

% --- Roof  parameters----------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------
A1=2.905392768;
A2=-3.204944487;
A3=1.320700525;
phi_max=0.45; % Kg/m2;
n=40; 
T=691200; 
dt=36;
dt1=3600/dt;
m=T/dt;
a=0.2;
Tini=299;  % K 
L=2257000;         % Latent heat of vaporization of water, J/Kg K
alpha = 0.35;      % absorptivity of insulation
sigma=5.67e-08;    % stefan boltzman constant
emissivity=0.94;   % emissivity of Insulation
lewis=0.83;           % lewis number
rhoda=1.2041;      % density of dry air
phiini=0.06; % initial water content

mFileErrorCode = 110    % After setting parameters


% --- Process Inputs ---------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------

sw1   = trnInputs(1);
RH1   = trnInputs(2);
Tout1 = trnInputs(3);
p1    = trnInputs(4);
v1    = trnInputs(5);
Troom = trnInputs(6);
alpharoom = trnInputs(7);
Qird  = trnInputs(8);
Pr    = trnInputs(9);

mFileErrorCode = 120    % After processing inputs


% --- First call of the simulation: initial time step (no iterations) --------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------
% (note that Matlab is initialized before this at the info(7) = -1 call, but the m-file is not called)

if ( (trnInfo(7) == 0) & (trnTime-trnStartTime < 1e-6) )  
    
    % This is the first call (Counter will be incremented later for this very first call)
    nCall = 0;

    % This is the first time step
    nStep = 1;
    
    % Initialize history of the variables for plotting at the end of the simulation
    nTimeSteps = (trnStopTime-trnStartTime)/trnTimeStep + 1;
    U=zeros(n+2,nTimeSteps);
    R=zeros(42,nTimeSteps);
    history.R       = zeros(42,nTimeSteps);
    history.CV      = zeros(nTimeSteps,1);
    history.CD      = zeros(nTimeSteps,1);
    history.Rnet    = zeros(nTimeSteps,1);
    history.TL      = zeros(nTimeSteps,1);
    history.Xsat    = zeros(nTimeSteps,1);
    history.Ttop    = zeros(nTimeSteps,1);
    history.Tbottom = zeros(nTimeSteps,1);
    history.U_tk    = zeros(n+2,1);
    history.phi     = zeros(nTimeSteps,1);
    dx=zeros(n,1);
    % creation of grid
for i=1:n
    dx(i)=0.2/n;
end
    U_tk=zeros(n+2,1);
    for i=1:n+2
        U_tk(i)=Tini;
    end
    U(:,1)=U_tk;
    phi_tk=phiini;
 % definition of thermo physical properties
    lambda=zeros(n,1);
    rho=zeros(n+2,1);
    cp=zeros(n+2,1);
    
   for i=1:n
        lambda(i)= 2.1;
        rho(i)=2400;
        cp(i)=1000;
   end
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
% complete definition of grid and Time
t1=linspace(0,T,nTimeSteps);
  x=zeros(n+2,1);
   x(1)=0.0;
   x(2)=dx(1)/2;
   for i=2:40
      x(i+1)=x(i)+0.5*(dx(i-1)+dx(i));
   end
   x(n+2)=a;
   
   D= zeros(n+2,2);
    
     Theta=zeros(2,1);
     F=zeros(n+2,1);
     
    mFileErrorCode = 130    % After initialization
    
end


% --- Very last call of the simulation (after the user clicks "OK"): Do nothing ----------------------------------------
% ----------------------------------------------------------------------------------------------------------------------

if ( trnInfo(8) == -1 )

    mFileErrorCode = 1000;
    
    % Draw a plot of Temperature profile along depth
    plot(x,history.U_tk(:,nTimeSteps));
    title('variation of temperature');
    ylabel('Temperature [K]');
    xlabel('depth [m]');
   
    mFileErrorCode = 0; % Tell TRNSYS that we reached the end of the m-file without errors
    return

end


% --- Post convergence calls: store values -----------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------

if (trnInfo(13) == 1)
    
    mFileErrorCode = 140;   % Beginning of a post-convergence call 
    history.U_tk(:,nStep)= U_tk_1;
    history.phi(nStep)=phi_tk;
 
    U_tk=U_tk_1;
    
     if phi_tk_1<= phi_max
            phi_tk=phi_tk_1;
        else
            phi_tk_1=phi_max;
            phi_tk=phi_tk_1;
     end
 
       
     
    
    mFileErrorCode = 0; % Tell TRNSYS that we reached the end of the m-file without errors
    return  % Do not update outputs at this call

end


% --- All iterative calls ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------

% --- If this is a first call in the time step, increment counter ---

if ( trnInfo(7) == 0 )
    nStep = nStep+1;
end

% --- Get TRNSYS Inputs ---

nI = trnInfo(3);     % For bookkeeping
nO = trnInfo(6);   % For bookkeeping

sw1   = trnInputs(1);
RH1   = trnInputs(2);
Tout1 = trnInputs(3);
p1    = trnInputs(4);
v1    = trnInputs(5);
Troom = trnInputs(6);
alpharoom = trnInputs(7);
Qird  = trnInputs(8);
Pr    = trnInputs(9);

mFileErrorCode = 150;   % After reading inputs 

% --- Main loop ---
       alphao=6.42+3.96*v1;
     
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
    D(1,1)= alphao;
    D(n+2,2)= alpharoom;
    Theta(1)=Tout1+273.15;
    Theta(2)=Troom+273.15;
    E=D*Theta; 
   % creation of main loop
        Pwss=6.116441*10^(7.591386*(U_tk(1)-273.15)/((U_tk(1)-273.15)+240.7263));
        Xss=0.622*Pwss/(p1-Pwss);
       Pwso=6.116441*10^(7.591386*Tout1/(Tout1+240.7263));
        Pwo=RH1*Pwso/100;
        Xo=0.622*Pwo/(p1-Pwo);
        f=Pwo*10^2/(13.6*9.81);
        epsky= 0.526+0.076*sqrt(f);
        Isky=epsky*sigma*(Tout1+273.15)^4;
        Ilw=emissivity*Isky-emissivity*sigma*(U_tk(1))^4;
        Rnet=alpha*sw1+emissivity*Isky-emissivity*sigma*(U_tk(1))^4;
        Cpm=1000+Xo*2100;
        rhoma=rhoda*(1+Xo)/(1+1.609*Xo);
        Kx=alphao/(lewis*rhoma*Cpm);
        phi_tk_1=phi_tk+dt*(Pr-((A1*(phi_tk/phi_max)^3+A2*(phi_tk/phi_max)^2+A3*(phi_tk/phi_max))*Kx*(Xss-Xo)));
        beta=A1*(phi_tk/phi_max)^3+A2*(phi_tk/phi_max)^2+A3*(phi_tk/phi_max);
        EV=beta*Kx*(Xss-Xo);
        F(1,1)=Rnet-L*EV;
        F(n+2,1) =Qird*1.8;
        Z=(M1*U_tk)+E+F;
        U_tk_1=A\Z;
        phi=phi_tk;
        CV=alphao*((Tout1+273.15)-U_tk(1));
        CD=c(1)*(U_tk(2)-U_tk(1));
        TL=alpharoom*(U_tk(n+2)-(Troom+273.15));
        Ttop=U_tk(1);
        Tbottom=U_tk(n+2);
% --- Set outputs ---

trnOutputs(1) = Ttop;
trnOutputs(2) = Tbottom;
trnOutputs(3) = CD;
trnOutputs(4) = CV;
trnOutputs(5) = TL;
trnOutputs(6) = Rnet;
trnOutputs(7) = Isky;
trnOutputs(8) = Ilw;
trnOutputs(9) = phi;
trnOutputs(10) = L*EV;
mFileErrorCode = 0; % Tell TRNSYS that we reached the end of the m-file without errors
return