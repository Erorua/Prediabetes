function GB_Bifurcation_Glucotoxicity


clc;


% ============================================================================================
% Description
% ============================================================================================


%%% Author: Aurore Woller

%%% Date: December 2022

%%% Uni: Weizmann institute of Science

%%% Description:Mathem. Model with carrying capacity C and glucotoxicity

%%% Simulates  Si, G and Beta vs time


%%% For Si decreasing exponantially with time as Si=Si0*exp(alpha*t)







% ============================================================================================
% Main
% ============================================================================================


%%% Kinetic param (in h and not in min)


u0=48;



Sip=1.5;

C=1.44;

a=7.85;




mu=0.05;


G1=23;

G0=4.8;

K=1e+05;



th=500;



slope=1e-03;


y0=1.5;


Si0=y0*exp(slope*th);



kine=[u0,C,Sip,a,G0,G1,K,mu,th,slope,y0,Si0];





%%%% Task:

tdyn(kine);




%====================================================================
% Time dynamics
%====================================================================

function tdyn(kine);








%%% kinetic parameters
u0=kine(1);
C= kine(2);
Sip=kine(3);
a=kine(4);
G0=kine(5);
G1=kine(6);
K=kine(7);
mu=kine(8);
th=kine(9);
slope=kine(10);
y0=kine(11);
Si0=kine(12);
% 


ratio=[0.01 0.1 1 10 100];


 


CC = {[0.05,0.39,0.16],[0.14,0.76,0.64],[0.00,1.00,1.00],[0.07,0.62,1.00],[0.00,0.00,1.00]} 




i=4;




%%%% Initial conditions:





g0=5;%- after bifurc


fg0=(g0^2)/((a^2)+(g0^2));



beta0=(u0/g0-C)*1/((fg0*y0/432));


%beta0=30


v=[g0 beta0];


%%%% Time:

trans=20000;
tend=20000;
tstep=0.01;




%%%% Transient Integration:

ttrans = [0:tstep:trans];
tspan = [0:tstep:tend];

option = odeset('RelTol', 1e-5);

if trans > 0 
    
 th=100100; 
 
 kine(9)=th;
    


[t x] = ode45(@dxdt,ttrans,v,option,kine); % transient integration: no HF, no jetlag
v=x(end,:);


end

 th=500; 
 
 kine(9)=th;
    


[t x] = ode45(@dxdt,tspan,v,option,kine);


%  Formula for change in SI

Sip=(y0).*heaviside(th-t)+heaviside(t-th).*((Si0)*exp(-slope*t));




figure(1)


plot(t/365,Sip,'.','MarkerSize',10,'color',CC{i})

xlabel('Time (days)')

ylabel('Insulin sensitivity (HOMA-IR-1)')


hold on;

set(gca,'FontName','Arial','FontSize',20);

xlim([0 12])

ylim([0.0 1.6])

pbaspect([1 1 1])








figure(2)


plot(t/365,(x(:,1)),'.','MarkerSize',10,'color',CC{i})




xlabel('Time (days)')

ylabel('Glucose (mmol/L)')

hold on;

set(gca,'FontName','Arial','FontSize',20);

pbaspect([1 1 1])


xlim([0 12])

ylim([0 15])


figure(3)


plot(t/365,(x(:,2)),'.','MarkerSize',10,'color',CC{i})


xlabel('Time (days)')

ylabel('Beta (\muU/ml.day-1)')

hold on;

set(gca,'FontName','Arial','FontSize',20);

pbaspect([1 1 1])


xlim([0 12])

ylim([0 2.5e+04])





% ============================================================================================
% dvdt
% ============================================================================================

function dv = dxdt(t,v,kine)


%%% variables

G=v(1);
beta=v(2);


%%% kinetic parameters
u0=kine(1);
C= kine(2);
Sip=kine(3);
a=kine(4);
G0=kine(5);
G1=kine(6);
K=kine(7);
mu=kine(8);
th=kine(9);
slope=kine(10);
y0=kine(11);
Si0=kine(12);





%  Formula for change in SI

Sip=(y0).*heaviside(th-t)+heaviside(t-th).*((Si0)*exp(-slope*t));



 

% Hill function

fg=(G^2)/((a^2)+(G^2));

% Kinetic parameters


G1=23;


gamma=432;

K=1e+05; % carrying capacity 


% %%% equations

% % 
dv = [
    u0-((fg*Sip*beta/gamma)+C)*G; % dg/dt
    beta*(mu)*((G/G0+G/G1)*(1-(beta/K)^1)-1-(G^2)/(G0*G1));    % dbeta/dt

] ; 



