function[G_EDO_Global,para_EDO,t,pars,pars_ci,Scale_ODE,Sc_ODE,Sm_ODE]=...
    ODE_optimi(T,Nt,dt,x_data,y_data,InitialMerozoites,r0,NbStage,beta,...
    Delta0,mStarM,mStarC,gamma_r,gamma_m,gamma_s,mu_m,mu_sd,mu_rm,mu_ms,...
    Lambda0,r)

%{
Understanding dynamics of P falciparum
gametocytes production: Insights from an age-structured model

Fitting and solving the ODE model
T= time of integration
Nt= number of points on [0,T]
dt= time step
IR parameters: Delta0; mStarM; mStarC;
mu_m=Decay rates of malaria parasites;
mu_x=natural death rates for uninfected RBCs;
mu_sd;mu_rm;mu_ms= duration of the RBC stages;
r=Merozoites multiplication factor;
mu_bar=parameter for the function mu(a);
Lambda0=Production rate of RBC;
gamma_r;gamma_m;gamma_s= RBCs preference coef;
        
        %}
    

Tspan=x_data;
function yu=myfun(param,Tspan)
    %{Here is the function to optimized: Model VS Data%}
mu_p=param(1);mu_g=param(2);
Sam=0;

%Cerating state variables (as vectors with length= length(data))  
Rr=x_data; Rm=x_data; Rs=x_data; 
P=x_data; m=x_data;para=m;
Sc=m;Sm=m;G=m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solving the ODE with an implicite volume scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initializing state variables
Rr(1)=Lambda0/mu_rm; Rm(1)=Lambda0/mu_ms; Rs(1)=Lambda0/mu_sd; m(1)=InitialMerozoites;
P(1)=0;
para(1)=100*P(1)/(P(1)+Rr(1)+Rm(1)+Rs(1));
G(1)=0;
%iterrating the system the ODE system with an implicte scheme to fit on data 
for n=1:length(x_data)-1
    dt1=24*(x_data(n+1)-x_data(n));%An adaptive time step in introduce to match data
    Dn0=floor(Delta0/dt1);
    ac=m(n)/(mStarC+m(n));Sc(n)=ac;
    if n<=Dn0
            am=0;
        elseif n>Dn0
            Sam=Sam+dt1*m(n-Dn0); 
            am=Sam/(mStarM+Sam);Sm(n)=am;
    end
    
    Rr(n+1)=(Rr(n)+dt1*Lambda0)/(1+dt1*(mu_rm+beta*gamma_r*m(n)));
    Rm(n+1)=(Rm(n)+dt1*mu_rm*Rr(n))/(1+dt1*(mu_ms+beta*gamma_m*m(n)));
    Rs(n+1)=(Rs(n)+dt1*mu_ms*Rm(n))/(1+dt1*(mu_sd+beta*gamma_s*m(n)));
    P(n+1)=(P(n)+dt1*beta*(gamma_r*Rr(n)+gamma_m*Rm(n)+gamma_s*Rs(n))*m(n))/(1+dt1*mu_p);
    m(n+1)=(m(n)+dt1*(r0*r*mu_p*P(n)))/(1+dt1*(mu_m+(ac+am)+beta*(gamma_r*Rr(n)+gamma_m*Rm(n)+gamma_s*Rs(n))));
       
    para(n+1)=100*P(n+1)/(P(n+1)+Rr(n+1)+Rm(n+1)+Rs(n+1));

    IIg=(1-r0)*r*mu_p*P(n);
    G(n+1)=(G(n)+dt1*IIg)/(1+dt1*mu_g);
end
yu=(max(y_data)/max(G))*G;%output variable to be optimized (VS the data)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nonlinear fit of the function 'myfun' defined above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter to estimate (mu_p,mu_g) is this order: 
%mu_p=death rate of pRBCs; mu_g=Proportion of sexual merozoites;
parguess =[1/48,1/500.10]; % nonlinear fits need initial guesses
lb=[1/48,0];ub=[1/48,1];% lower and upper bonds
[pars,resnorm,residual,exitflag,output_fit,lambda_fit,jacobian]=...
    lsqcurvefit(@myfun,parguess,x_data,y_data,lb,ub);
alpha = 0.005; % this is for 95% confidence intervals
pars_ci = nlparci(pars,residual,'jacobian',jacobian,'alpha',alpha);
%estimated parameters are given by the variable 'pars'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%we now solve the ODE systems over time T, for a given k=NbStage pRBCs stage, with estimated parameters above 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global EDO solution (with the implicit volume scheme)
param=pars;%here are parameters estimated above
mu_p=param(1);mu_g=param(2);

%creating state variables
Sam=0;
Rr=zeros(1,Nt+1); Rm=zeros(1,Nt+1); Rs=zeros(1,Nt+1); 
P=zeros(1,Nt+1); m=zeros(1,Nt+1); para=m;
PP=zeros(NbStage,Nt+1);
Sc=m;Sm=m;G=m;

%initila conditions
Rr(1)=Lambda0/mu_rm; Rm(1)=Lambda0/mu_ms; Rs(1)=Lambda0/mu_sd; m(1)=InitialMerozoites;
P(1)=0;
for j_stage=1:NbStage, PP(j_stage,1)=0; end
para(1)=100*P(1)/(P(1)+Rr(1)+Rm(1)+Rs(1));
G(1)=0;
%integrating the ODE system with K stages
for n=1:Nt
    Dn0=floor(Delta0/dt);
    ac=m(n)/(mStarC+m(n));Sc(n)=ac;
    if n<=Dn0
            am=0;
        elseif n>Dn0
            Sam=Sam+dt*m(n-Dn0);
            am=Sam/(mStarM+Sam);Sm(n)=am;
    end
    
    Rr(n+1)=(Rr(n)+dt*Lambda0)/(1+dt*(mu_rm+beta*gamma_r*m(n)));
    Rm(n+1)=(Rm(n)+dt*mu_rm*Rr(n))/(1+dt*(mu_ms+beta*gamma_m*m(n)));
    Rs(n+1)=(Rs(n)+dt*mu_ms*Rm(n))/(1+dt*(mu_sd+beta*gamma_s*m(n)));
    PP(1,n+1)=(PP(1,n)+dt*beta*(gamma_r*Rr(n)+gamma_m*Rm(n)+gamma_s*Rs(n))*m(n))/(1+dt*NbStage*mu_p);
    for j_stage=2:NbStage
        PP(j_stage,n+1)=(PP(j_stage,n)+dt*NbStage*mu_p*PP(j_stage-1,n))/(1+dt*NbStage*mu_p);
    end
    P(n+1)=PP(NbStage,n+1);
    m(n+1)=(m(n)+dt*(r0*r*NbStage*mu_p*P(n)))/(1+dt*(mu_m+(ac+am)+beta*(gamma_r*Rr(n)+gamma_m*Rm(n)+gamma_s*Rs(n))));
    
    
para(n+1)=100*sum(PP(:,n+1))/(sum(PP(:,n+1))+Rr(n+1)+Rm(n+1)+Rs(n+1));

IIg=(1-r0)*r*NbStage*mu_p*P(n);
G(n+1)=(G(n)+dt*IIg)/(1+dt*mu_g);
end

%some outputs
Scale_ODE=(max(y_data)/max(G));
G_EDO_Global=Scale_ODE*G;
para_EDO=para;
Sc_ODE=Sc;
Sm_ODE=Sm;
%%End global solution EDO


t=0:dt:T;t=t/24;
end
