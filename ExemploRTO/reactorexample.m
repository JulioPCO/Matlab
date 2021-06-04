clc 
clear all
close all

addpath ('C:\Users\julio\OneDrive\Área de Trabalho\Project Norway\Casadi')
import casadi.*

%% Parameters
%reaction rates
% k1 = 0.75*10^4;  %[m3/(kmol min)^2]
% k2 = 1.5*10^4;  %[m3/(kmol min)]
V = 50; %[m^3]
rhoab = 7; %[kmol/m3] 
rhoc = 5.5; %[kmol/m3]
Fabmol = 1000; %[kmol/min]

%% Differential states
%symbolic declaration
%mole fractions 
% reactor = [Ca, Cc, Cd, Ce]
x_var = MX.sym('x_var',4,1); %[kmol/m3]

%% Input
%xa(nt), Fpurge, Fc (volumetric)
u_var = MX.sym('u_var',3,1); %[m3/min]

%time transformation: CASADI integrates always from 0 to 1 and the USER does the time
%scaling with T.
T = MX.sym('T',1); %[min]

%Symbolic reaction rate
k1 = MX.sym('k1',1); %[m3/(kmol min)^2]
k2 = MX.sym('k2',1); %[m3/(kmol min)]

%inputs and parameters
p_var = [k1; k2; T];

%% System equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Reactor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reaction array
%A + C -> D
r1 = k1*10^4*x_var(1)*x_var(2);
%2C -> E
r2 = k2*10^4*x_var(2)^2;

%mass balance
Fab = Fabmol/rhoab - u_var(2);%volumetric, density is constant
Fout = Fab + u_var(3);%Fout = Fd + Ffr %[m3/min]

%component balances
dCadt = (Fab*(u_var(1)*rhoab) - Fout*x_var(1))/V - r1; %C_A
dCcdt = (u_var(3)*rhoc - Fout*x_var(2))/V - r1 - 2*r2;%C_C pure C is fed to the system
dCddt = -Fout*x_var(3)/V + r1;%C_D
dCedt = -Fout*x_var(4)/V + r2;%C_E

%end modeling

%% give parameter values
% Form the ODE system
diff = [dCadt;dCcdt;dCddt;dCedt];


%% Casadi commands

   % ================================================
   %       Integrator
   % ================================================

%declaring function in standard DAE form (scaled time)
ode = struct('x',x_var,'p',vertcat(u_var,p_var),'ode',T*diff);
      
%calling the integrator, the necessary inputs are: label; integrator; function with IO scheme of a DAE (formalized); struct (options)
F = integrator('F','cvodes',ode);


%initial guess
dx0 = [2.2; 1; 0; 0];
u0 = [0.91;530/rhoab;400/rhoc];
k10 = 0.75;
k20 = 1.5;
Time = 0.02; %[min]
par0=[k10; k20; Time];
keepx=[dx0];
keepISS = [];
keepI = [];
keepestimate = [];
keepoptimum = [];

m=1;

for i = 1:1000

if i == 500
   u0 = [0.75;530/7;400/5.5];
end

Results = F('x0',dx0, 'p', [u0;par0]);
dx0 = [full(Results.xf(1)); full(Results.xf(2));full(Results.xf(3:4))] + 0.003*randn(4,1);
keepx = [keepx dx0];

if i == 50*m
   [ISS] = FtestSS(keepx(:,50*(m-1)+1:50*m)) ;
   keepISS = [keepISS ISS];
   
   [I,estimate] = Estimation(keepx(:,50*(m-1)+1:50*m),ISS,u0);
   keepI = [keepI;I];
   keepestimate = [keepestimate; estimate];
   
   [optimum] = Optimization(estimate,u0);
   keepoptimum = [keepoptimum optimum];
   
   if keepoptimum(1,m)>0
       u0(2:3)=optimum;
   end
   m=m+1
end

end

figure(8)
plot([(1:i)]/100, keepx(:,1:i))
title('Concentration x time')
legend('Ca','Cc','Cd','Ce')
xlabel('time [min]')
ylabel('Concentration [kmol]')

            