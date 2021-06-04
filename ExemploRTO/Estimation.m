function [I,estimate] = Estimation(keepx, ISS,uinit)
addpath ('C:\Users\julio\OneDrive\�rea de Trabalho\Project Norway\Casadi')
import casadi.*

%estimate initialization
estimate = -1;
%global for steady state
global ii
ii=1;

%Finding 10 steady state point consecutive
i=1;
count=0;
I=0;

while i <= size(ISS,1)
    if ISS(i,1) == 1
        count=count+1;
    else
        count=0;
    end
    if count==15
        ii=i-9;
        i=size(ISS,1);
        I=1;
    end
    i=i+1;
end

if count < 15
   return  
end

xSS = [];
for j=1:4
    xSS = [xSS; mean(keepx(j,ii:ii+9))];
end

%% Model
%% Parameters
%reaction rates
% k1 = 0.75*10^4;  %[m3/(kmol min)^2]
% k2 = 1.5*10^4;  %[m3/(kmol min)]
V = 50; %[m3]
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
rhoab = MX.sym('rhoab',1); %[m3]

%inputs and parameters
p_var = [k1; k2; rhoab;T];

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
%% Estimation Problem
% ===================================
%     Declaring 
% ===================================
%Vector and symbolic values for x and p
xp={};
xp0=[];
xplb=[]; %lower bound constraints for x and p
xpub=[];

%Vector and symbolic value for constrainst referring the 
%differential equation
g={};
gub=[];
glb=[];

%initial guess
dx0 = [2.2; 1; 0; 0];
u0 = uinit;
k10 = 0.75;
k20 = 1.5;
rhoab0 = 7; %[kmol/m3] 
Time = 0.02; %[min]
par0=[k10; k20; rhoab0;Time];

%Constraints 
lbx = [0;0;0;0];
ubx = [5;5;5;5];
lbp = [k10;k20;5;Time];
ubp = [k10;k20;10;Time];

%Assigning to vectors
xp={xp{:},x_var,u_var,p_var};
xp0=[xp0;dx0;u0;par0];
xplb=[xplb; lbx;uinit;lbp];
xpub=[xpub;ubx;uinit;ubp];

g={g{:}, vertcat(diff)};
glb=[glb;zeros(4,1)];
gub=[gub;zeros(4,1)];


%objective function - total modifier
J =  (x_var(1) - xSS(1))'*(x_var(1) - xSS(1)) +(x_var(3) - xSS(3))'*(x_var(3) - xSS(3))+(x_var(4) - xSS(4))'*(x_var(4) - xSS(4));

%Solver Options
options = struct;
options.ipopt.acceptable_tol = 1e-8;


%Nlp structure and solver
nlp=struct('x',vertcat(xp{:}),'f',J,'g',vertcat(g{:}));
solver=nlpsol('solver','ipopt',nlp,options);
sol=solver('x0',xp0,'lbx',xplb,'ubx',xpub,'lbg',glb,'ubg',gub);

if solver.stats.success ~=1
    e=1;
    return
end
% clc 
%clean ipopt output

% Extract Solution
objfun = full(sol.f)
estimate = full(sol.x(10))
end