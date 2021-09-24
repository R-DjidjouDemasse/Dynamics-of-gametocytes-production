function[G_EDO_Global,para_EDO,t,pars,pars_ci,Scale_ODE,NormDiffEDO,Sc_ODE,Sm_ODE]=...
    ODE_optimi(x_data,y_data,InitialMerozoites,r0,NbStage,beta,...
    Delta0,mStarM,mStarC,gamma_r,gamma_m,gamma_s,mu_m,mu_sd,mu_rm,mu_ms,...
    Lambda0,r)

%{
Understanding dynamics of P falciparum
gametocytes production: Insights from an age-structured model

Fitting and solving the ODE model

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
mu_p=param(1);mu_g=param(2);
Sam=0;

    

    
Rr=x_data; Rm=x_data; Rs=x_data; 
P=x_data; m=x_data;para=m;
Sc=m;Sm=m;G=m;

%solving the ODE with an implicite volume scheme
Rr(1)=Lambda0/mu_rm; Rm(1)=Lambda0/mu_ms; Rs(1)=Lambda0/mu_sd; m(1)=InitialMerozoites;
P(1)=0;
para(1)=100*P(1)/(P(1)+Rr(1)+Rm(1)+Rs(1));
G(1)=0;
for n=1:length(x_data)-1
    dt=24*(x_data(n+1)-x_data(n));Dn0=floor(Delta0/dt);
    ac=m(n)/(mStarC+m(n));Sc(n)=ac;
    if n<=Dn0
            am=0;
        elseif n>Dn0
            Sam=Sam+dt*m(n-Dn0);
            %Nm=floor((n*dt-8*24)/dt);
            %am=dt*trapz(m(1:Nm)); 
            am=Sam/(mStarM+Sam);Sm(n)=am;
    end
    
    Rr(n+1)=(Rr(n)+dt*Lambda0)/(1+dt*(mu_rm+beta*gamma_r*m(n)));
    Rm(n+1)=(Rm(n)+dt*mu_rm*Rr(n))/(1+dt*(mu_ms+beta*gamma_m*m(n)));
    Rs(n+1)=(Rs(n)+dt*mu_ms*Rm(n))/(1+dt*(mu_sd+beta*gamma_s*m(n)));
    P(n+1)=(P(n)+dt*beta*(gamma_r*Rr(n)+gamma_m*Rm(n)+gamma_s*Rs(n))*m(n))/(1+dt*mu_p);
    m(n+1)=(m(n)+dt*(r0*r*mu_p*P(n)))/(1+dt*(mu_m+(ac+am)+beta*(gamma_r*Rr(n)+gamma_m*Rm(n)+gamma_s*Rs(n))));
    
    
para(n+1)=100*P(n+1)/(P(n+1)+Rr(n+1)+Rm(n+1)+Rs(n+1));

IIg=(1-r0)*r*mu_p*P(n);
%G(n+1)=G(n)+dt*(IIg-mu_g*G(n));
G(n+1)=(G(n)+dt*IIg)/(1+dt*mu_g);
end
yu=(max(y_data)/max(G))*G;
end


parguess =[1/48,1/500.10]; % nonlinear fits need initial guesses
%x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
lb=[1/48,0];ub=[1/48,1];
[pars,resnorm,residual,exitflag,output_fit,lambda_fit,jacobian]=...
    lsqcurvefit(@myfun,parguess,x_data,y_data,lb,ub);
alpha = 0.005; % this is for 95% confidence intervals
pars_ci = nlparci(pars,residual,'jacobian',jacobian,'alpha',alpha);

%%%%%GLOBAL EDO SOLUTION (with an implicit volume scheme)
pars(2)=1/(1/pars(2)+x_data(7)*24);
param=pars;
T=40*24; Nt=T/20;
dt=T/Nt;

mu_p=param(1);mu_g=param(2);
 
Sam=0;
Rr=zeros(1,Nt+1); Rm=zeros(1,Nt+1); Rs=zeros(1,Nt+1); 
P=zeros(1,Nt+1); m=zeros(1,Nt+1); para=m;
PP=zeros(NbStage,Nt+1);
Sc=m;Sm=m;G=m;
 
Rr(1)=Lambda0/mu_rm; Rm(1)=Lambda0/mu_ms; Rs(1)=Lambda0/mu_sd; m(1)=InitialMerozoites;
P(1)=0;
for j_stage=1:NbStage, PP(j_stage,1)=0; end
para(1)=100*P(1)/(P(1)+Rr(1)+Rm(1)+Rs(1));
G(1)=0;
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
    m(n+1)=(m(n)+dt*(r0*r*mu_p*P(n)))/(1+dt*(mu_m+(ac+am)+beta*(gamma_r*Rr(n)+gamma_m*Rm(n)+gamma_s*Rs(n))));
    
    
para(n+1)=100*sum(PP(:,n+1))/(sum(PP(:,n+1))+Rr(n+1)+Rm(n+1)+Rs(n+1));

IIg=(1-r0)*r*mu_p*P(n);
G(n+1)=(G(n)+dt*IIg)/(1+dt*mu_g);
end
Scale_ODE=(max(y_data)/max(G));
G_EDO_Global=Scale_ODE*G;
para_EDO=para;
Sc_ODE=Sc;
Sm_ODE=Sm;

yu1=2:length(x_data);
for n1=2:length(x_data)
    t1=24*x_data(n1); nt1=floor(t1/dt);
    yu1(n1)=(G_EDO_Global(nt1)+(t1-nt1)*(G_EDO_Global(nt1)-G_EDO_Global(nt1+1))/dt);
end
yu1=yu1-y_data'; NormDiffEDO=norm(yu1,2);
%%END GLOBAL SLOTUION EDO


t=0:dt:T;t=t/24;
end
