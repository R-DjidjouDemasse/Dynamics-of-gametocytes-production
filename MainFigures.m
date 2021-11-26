function[Name,mu_g,AlphaG_PDE,AlphaG_ODE,InitialMerozoites,...
    mu_p_EDO,mu_g_EDO]=MainFigures

%{
Understanding dynamics of P falciparum
gametocytes production: Insights from an age-structured model
    
Solving the ODE and PDE models

Required the function 'ODE_optimi'
    function[G_EDO_Global,para_EDO,t,pars,pars_ci,Scale_ODE,NormDiffEDO,Sc_ODE,Sm_ODE]=...
    ODE_optimi(T,Nt,dt,x_data,y_data,InitialMerozoites,r0,NbStage,beta,...
    Delta0,mStarM,mStarC,gamma_r,gamma_m,gamma_s,mu_m,mu_sd,mu_rm,mu_ms,...
    Lambda0,r)

IR parameters: Delta0; mStarM; mStarC;
mu_m=Decay rates of malaria parasites;
mu_x=natural death rates for uninfected RBCs;
mu_sd;mu_rm;mu_ms= duration of the RBC stages;
r=Merozoites multiplication factor;
mu_bar=parameter for the function mu(a);
Lambda0=Production rate of RBC;
gamma_r;gamma_m;gamma_s= RBCs preference coef;
        
        %}

%Here are some fixed parameters
Delta0=16*24; mStarM=20.4*10^6; mStarC=2755*10^3;
mu_m=48/24;mu_x=0.00833/24;
mu_sd=1/(48);mu_rm=1/(36);mu_ms=1/(116.5*24);
r=16;
r0=0.95; 
mu_bar=10;
Lambda0=1.73*10^6;

T=40*24;%Time of integration
Nt=10*T;%Number of points on [0,T]
dt=T/Nt;%Time step. Here, Nt is set such that the time step dt=0.1. 
%Small dt is necessary for the stability of the implicit/explicit scheme

a_max=2.47*24;%max age of pRBCS
Na=590;%number of points on [0,a_max]
da=a_max/Na;%age step

%RBCs preference for P. Falciparum
gamma_r=1; gamma_m=1; gamma_s=1;

%a given data set
Data=xlsread('Patient_Data.xlsx');
seques_time=48;
InitialMerozoites=10^7;
Name='Patient name';

%vectors for the number K of stages for the EDO model
VectNbStage=[1,40,100,150];

beta=(6.2734*10^-9)/24;%Infection rate of uRBC

%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving the ODE systems with one stage (k=1)
%%%%%%%%%%%%%%%%%%%%%%%%%%
x_data=Data(:,1); y_data=Data(:,2); 
NbStage=VectNbStage(1);%the number of stage (k=1)
%Fitting and solving the ODE model with the function 'ODE_optimi'
[G_EDO_Global,para_EDO,t_EDO,pars,pars_ci,Scale_ODE,Sc_ODE,Sm_ODE]=...
    ODE_optimi(T,Nt,dt,x_data,y_data,InitialMerozoites,r0,NbStage,beta,...
    Delta0,mStarM,mStarC,gamma_r,gamma_m,gamma_s,mu_m,mu_sd,mu_rm,mu_ms,...
    Lambda0,r);
%here are estimated parameters
mu_p_EDO=pars(1);
mu_g_EDO=pars(2);
%NOTE: above estimated parameters for the ODE are used for the PDE system
mu_g=mu_g_EDO;
AlphaG_ODE = (1-r0)*Scale_ODE;
r0=1-AlphaG_ODE/Scale_ODE;


%For figure's legend etc ....
PlotStyle={'--k',':k','-.k','--r'};
TextLegendGameto={};
TextLegendGameto{1}='PDE prediction';
for iv=1:length(VectNbStage)
    TextLegendGameto{iv+1}=['ODE prediction $(K=',num2str(VectNbStage(iv)),')$'];
end
TextLegendGameto{iv+2}='Observed data';
TextLegendParaEDO={};
for iv=1:length(VectNbStage)
    TextLegendParaEDO{iv}=['Parasitemia$(K=',num2str(VectNbStage(iv)),')$'];
    TextLegendParaEDO{iv+length(VectNbStage)}=['Gametocytemia$(K=',num2str(VectNbStage(iv)),')$'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving the PDE system with an Euler explicit method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dn0=floor(Delta0/dt);Sam=0;

%RBC'srupture function
function[y]=mu(a)
    if a<seques_time
        y=0;
    else, y=mu_bar;
    end
end

%State system for th PDE model
function[y]=state(Lambda,Rr,Rm,Rs,m,P)


    Mat1=[0,0,0,0,0;1,-mu_rm,0,0,0;0,mu_rm,-mu_ms,0,0;0,0,mu_ms,-mu_sd,0;0,0,0,0,-mu_m-(ac+am)];
     Mat2=zeros(Na,Na);
     for j=2:Na
         Mat2(j,j)=-(1/da+mu(j*da)+mu_x);
         Mat2(j,j-1)=1/da;
     end
     Mat2(1,1)=-(1/da+mu(da)+mu_x);
     Mat=zeros(Na+5,Na+5);
     Mat(1:5,1:5)=Mat1; Mat(6:Na+5,6:Na+5)=Mat2;
    phi=Lambda0-Lambda+mu_sd*Rs+beta*(gamma_r*Rr+gamma_m*Rm+gamma_s*Rs)*m;
    FF=zeros(5+Na,1);
    FF(1)=0; FF(2)=-beta*gamma_r*Rr*m; FF(3)=-beta*gamma_m*Rm*m; FF(4)=-beta*gamma_s*Rs*m;
     P1=[beta*(gamma_r*Rr+gamma_m*Rm+gamma_s*Rs)*m;P];
    %yy1=zeros(1,Na);
    s1=r*mu(0)*beta*(gamma_r*Rr+gamma_m*Rm+gamma_s*Rs)*m...
  +r*mu(da*Na)*P1(Na+1)+4*r*mu((Na-1)*da)*P1(Na);
    for j=1:Na/2
        s1=s1+2*r*mu(2*j*da)*P1(2*j+1)+4*r*mu((2*j-1)*da)*P1(2*j);
    end
    FF(5)=r0*s1*da/3-beta*(gamma_r*Rr+gamma_m*Rm+gamma_s*Rs)*m;
     FF(6)=1/da*beta*(gamma_r*Rr+gamma_m*Rm+gamma_s*Rs)*m;
    y=Mat*[Lambda;Rr;Rm;Rs;m;P]+FF;
end

%cerating state variables
Lambda=zeros(1,Nt+1); Rr=zeros(1,Nt+1); Rm=zeros(1,Nt+1); Rs=zeros(1,Nt+1); para1=Rs;
P=zeros(Na+1,Nt+1); m=zeros(1,Nt+1);para=m; ring=zeros(1,Nt+1);tropho=zeros(1,Nt+1);
shitz=zeros(1,Nt+1);Sc=m;Sm=m;G=m;

%Initalizing state variables
Lambda(1)=1.73*10^6; 
Rr(1)=Lambda0/mu_rm; Rm(1)=Lambda0/mu_ms; Rs(1)=Lambda0/mu_sd; m(1)=InitialMerozoites;
P(1,1)=beta*(gamma_r*Rr(1)+gamma_m*Rm(1)+gamma_s*Rs(1))*m(1);
s2=P(1,1)+P(Na+1,1)+4*P(Na,1);
for l=1:Na/2
    s2=s2+2*P(2*l+1,1)+4*P(2*l,1);
end
para(1)=s2*da/3; para1(1)=para(1);
para(1)=100*para(1)/(para(1)+Rr(1)+Rm(1)+Rs(1));
G(1)=0;

%Solving the PDE model, defined by the function 'state' above, using an
%Euler scheme
for n=1:Nt
    ac=m(n)/(mStarC+m(n));Sc(n)=ac;
    if n<=Dn0
            am=0;
        elseif n>Dn0
            Sam=Sam+dt*m(n-Dn0);
            am=Sam/(mStarM+Sam);Sm(n)=am;
    end
        
    k1=state(Lambda(n),Rr(n),Rm(n),Rs(n),m(n),P(2:Na+1,n));
    k2=state(Lambda(n)+dt/2*k1(1),Rr(n)+dt/2*k1(2),Rm(n)+dt/2*k1(3),Rs(n)+dt/2*k1(4),m(n)+dt/2*k1(5),P(2:Na+1,n)+dt/2*k1(6:Na+5));
    k3=state(Lambda(n)+dt/2*k2(1),Rr(n)+dt/2*k2(2),Rm(n)+dt/2*k2(3),Rs(n)+dt/2*k2(4),m(n)+dt/2*k2(5),P(2:Na+1,n)+dt/2*k2(6:Na+5));
    k4=state(Lambda(n)+dt*k3(1),Rr(n)+dt*k3(2),Rm(n)+dt*k3(3),Rs(n)+dt*k3(4),m(n)+dt*k3(5),P(2:Na+1,n)+dt*k3(6:Na+5));
    x1=[Lambda(n);Rr(n);Rm(n);Rs(n);m(n);P(2:Na+1,n)]+dt/6*(k1+2*k2+2*k3+k4);
    Lambda(n+1)=x1(1);Rr(n+1)=x1(2);Rm(n+1)=x1(3);Rs(n+1)= x1(4);
    m(n+1)=x1(5); P(2:Na+1,n+1)=x1(6:Na+5);
    
    P(1,n+1)=beta*(gamma_r*Rr(n+1)+gamma_m*Rm(n+1)+gamma_s*Rs(n+1))*m(n+1);
    s2=P(1,n+1)+P(Na+1,n+1)+4*P(Na,n+1);
    for l=1:Na/2
        s2=s2+2*P(2*l+1,n+1)+4*P(2*l,n+1);
    end
para(n+1)=s2*da/3;para1(n+1)=para(n+1);
para(n+1)=100*para(n+1)/(para(n+1)+Rr(n+1)+Rm(n+1)+Rs(n+1));

Ig=1:Na+1;
for J1=0:Na,Ig(J1+1)=r*mu(J1*da)*P(J1+1,n);end
IIg=(1-r0)*(da*trapz(Ig));
G(n+1)=(G(n)+dt*IIg)/(1+dt*mu_g);
end

t=0:dt:T; t=t/24;
%End solving the PDE model

%Defining the variable 'G_retard=G(t-seques_time)'
[Gmax11,IGmax11]=max(G);
time_delay=seques_time/dt;
G_retard=(time_delay+1):Nt+1;
for n2=(time_delay+1):Nt+1, G_retard(n2-time_delay)=G(n2);end



%Ploting Parasitemia and Gametocytemia
figure,
axes ('fontsize',14)
hold on
yyaxis left
plot(t,(max(Data(:,2))/max(G))*10^3*para,'LineWidth',2)
yyaxis right
plot(t,(max(Data(:,2))/max(G))*G,'--','LineWidth',2)
xlabel('Time (days)','fontsize',18,'Interpreter','latex'),
hold off
title(Name,'fontsize',18)
legend('Parasitemia (\%)','Gametocytemia per $\mu{L}$','fontsize',14,...
    'Interpreter','latex','location','southeast')

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ploting Gametocytemia for the PDE and ODE model in the same figures for
%different pRBC stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scale_PDE=(max(Data(:,2))/max(G));
figure,
axes ('fontsize',14)
hold on
%We first plot the PDE model's gameto
plot(t,Scale_PDE*G,'b','LineWidth',2),
%We now plot the ODE model's gameto (for k=1 pRBC stage)
plot(t_EDO,G_EDO_Global,PlotStyle{1},'LineWidth',2),
%Next, we solve and plot the ODE model' gameto for other values k (pRBC stage)
for IdStage=2:length(VectNbStage)
    NbStage=VectNbStage(IdStage);%number of pRBC stage
    %Solving the ODE model with the function 'ODE_optimi'
    [G_EDO_Global_0,para_EDO_0,t_EDO_0]=...
        ODE_optimi(T,Nt,dt,x_data,y_data,InitialMerozoites,r0,NbStage,beta,...
        Delta0,mStarM,mStarC,gamma_r,gamma_m,gamma_s,mu_m,mu_sd,mu_rm,mu_ms,...
        Lambda0,r);
    plot(t_EDO_0,G_EDO_Global_0,PlotStyle{IdStage},'LineWidth',2), 
end
plot(Data(:,1),Data(:,2),'o','LineWidth',2),
hold off
title(Name,'fontsize',18)
xlabel('Time (days)','fontsize',18,'Interpreter','latex'), 
ylabel('Mature Gametocyte per $\mu{L}$','fontsize',18,'Interpreter','latex')
legend(TextLegendGameto,...
    'fontsize',14,'location','southeast','Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ploting parasitemia for the PDE and ODE model in the same figure for
%different pRBC stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
axes ('fontsize',14)
hold on
%We first plot the PDE model's parasitemia
plot(t,(max(Data(:,2))/max(G))*10^3*para,'b','LineWidth',2),
%We now plot the ODE model's parasitemia (for k=1 pRBC stage)
plot(t_EDO,(max(Data(:,2))/max(G))*10^3*para_EDO,PlotStyle{1},'LineWidth',2),
%Next, we solve and plot the ODE model' parasitemia for other values k (pRBC stage)
for IdStage=2:length(VectNbStage)
    NbStage=VectNbStage(IdStage);%number of pRBC stage
    %Solving the ODE model with the function 'ODE_optimi'
    [G_EDO_Global_0,para_EDO_0,t_EDO_0]=...
        ODE_optimi(T,Nt,dt,x_data,y_data,InitialMerozoites,r0,NbStage,beta,...
        Delta0,mStarM,mStarC,gamma_r,gamma_m,gamma_s,mu_m,mu_sd,mu_rm,mu_ms,...
        Lambda0,r);
    plot(t_EDO_0,(max(Data(:,2))/max(G))*10^3*para_EDO_0,PlotStyle{IdStage},'LineWidth',2), 
end
hold off
title(Name,'fontsize',18)
xlabel('Time (days)','fontsize',18,'Interpreter','latex'), 
ylabel('Parasitemia (\%)','fontsize',18,'Interpreter','latex')
legend(TextLegendGameto,...
    'fontsize',14,'location','southeast','Interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for the link between Gameto and Parasitemia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first part of the estimate (i.e., for t<T0 in eq. 12 of the manuscript)
eG=log10((max(Data(:,2))/max(G))*G_retard(2:IGmax11)); 
eP=log10(para(2:IGmax11)); 
t_data=eP;
%Function to be fitted
function output=myfunc1(parameter,t_data)
    aaop(1)=parameter(1);
    aaop(2)=parameter(2);
    output=aaop(2)+ t_data*aaop(1);
end

opts = statset('nlinfit');
opts.RobustWgtFun = 'logistic';
pol1=polyfit(eP,eG,1);
parguess =[pol1(1),pol1(2)];
[aaop1, resid, J] = nlinfit(t_data,eG,@myfunc1,parguess,opts);
alpha = 0.05; % this is for 95% confidence intervals
pars_ci1 = nlparci(aaop1,resid,'jacobian',J,'alpha',alpha);
pars_ci1=pars_ci1-[aaop1(1) aaop1(1);aaop1(2) aaop1(2)];

apG=myfunc1(aaop1,t_data);
Rsq1 = 1 - sum((eG - apG).^2)/...
    sum((eG - mean(eG)).^2);

%second part of the estimate (i.e., for t>T0 in eq. 12 of the manuscript)
eG2=log10((max(Data(:,2))/max(G))*G(IGmax11:IGmax11+2250)); 
nnl2=length(IGmax11:IGmax11+2250);
eP2=-log10(para(IGmax11:IGmax11+2250)); 
t_data=eP2;
%Function to be fitted
function output=myfunc2(parameter,t_data)
    aaop(1)=parameter(1);
    aaop(2)=parameter(2);
    output=aaop(1)*t_data+aaop(2);
end
 pol=polyfit(eP2,eG2,1);
opts = statset('nlinfit');
opts.RobustWgtFun = 'logistic';
parguess2 =[pol(1),pol(2)];
[aaop2, resid, J] = nlinfit(t_data,eG2,@myfunc2,parguess2,opts);
alpha = 0.05; % this is for 95% confidence intervals
aaop2;
pars_ci2 = nlparci(aaop2,resid,'jacobian',J,'alpha',alpha);
pars_ci2=pars_ci2-[aaop2(1) aaop2(1);aaop2(2) aaop2(2)];
apG2=myfunc2(aaop2,t_data);
Rsq2 = 1 - sum((eG2(1:nnl2) - apG2(1:nnl2)).^2)/...
    sum(eG2(1:nnl2) - mean(eG2(1:nnl2)).^2);


figure,
axes ('fontsize',14)
hold on
plot(eP,eG,'b','LineWidth',2)%Ploting Parasitemia VS Gametocyte
plot(eP,apG,'--r','LineWidth',2)%Using an estimated formula
plot(-eP2,eG2,'b','LineWidth',2)%Ploting Parasitemia VS Gametocyte
plot(-eP2,apG2,'--r','LineWidth',2)%Using an estimated formula
hold off
title(Name,'fontsize',18)
ylabel('$Log_{10}$ Mature Gametocyte per $\mu{L}$','fontsize',18,'Interpreter','latex'), 
xlabel('$Log_{10}$ Parasitemia','fontsize',18,'Interpreter','latex')
legend('Ploting Parasitemia VS Gametocyte','Using an estimated formula',...
    'fontsize',14,'location','best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Printing some output on the sceen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AlphaG_PDE = (1-r0)*Scale_PDE;
fprintf(':::::::::::::::::::::::::::::::::::::::::%s.\n', ' '); 
fprintf('ID = %s.\n',Name); 
fprintf('mu_G = %d.\n',mu_g);
fprintf('AlphaG_PDE = %d.\n',AlphaG_PDE);
fprintf('InitialMerozoites = %d.\n',InitialMerozoites);


end
