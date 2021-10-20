function [Resqinit,Resqinit_e0,Resqinit_eps,OLSnew,OLSnew_e0,OLSnew_eps,veta1,veta,veps,b1,bL,b1_e0,bL_e0,b1_eps,bL_eps]=init()

global N K1 K2 K3 K3t K4 K5 K6 Ntau Vectau tau ...
       Y AGE meanAGE stdAGE meanY stdY ENT ...
       MatAGE1 AGE_t MatAGE_tot Y_t Y_lag Y_tot...
       Y1 AGE_tot Time_t K2t D_t T

 sample=find(D_t);
Time=[];
for tt=2:T
 Time=[Time tt*ones(N,1)]; 
end
Time_t=Time(:);
Time_t=Time_t(sample);

% Initial conditions: eta given eta_{t-1} and age  
Resqinit=zeros((K1+1)*(K2+1)*(K2t+1),Ntau);
MatYlag=[];
for kk1=0:K1
    for kk2=0:K2
        for kk2t=0:K2t
            MatYlag=[MatYlag hermite(kk1,(Y_lag+randn(size(Y_t,1),1)-meanY)/stdY).*hermite(kk2,(AGE_t-meanAGE)/stdAGE).*hermite(kk2t,Time_t)];
        end
    end
end

Ydep=Y_t+randn(size(Y_t,1),1);
for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(MatYlag,Ydep,tau);
    Resqinit(:,jtau)=beta1;
end
OLSnew=pinv(MatYlag)*Ydep;
veta=var(Ydep-MatYlag*OLSnew);

Vect1=Ydep-MatYlag*Resqinit(:,1);
Vect2=Ydep-MatYlag*Resqinit(:,Ntau);
b1=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
bL=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));


% Initial conditions: eta1 given age1 and period
Resqinit_e0=zeros((K4+1)*(K5+1)*(K6+1),Ntau);
Ydep = Y1 + 1*randn(N,1);
for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(MatAGE1,Ydep,tau);
    Resqinit_e0(:,jtau)=beta1;
end
OLSnew_e0=pinv(MatAGE1)*(Ydep);
veta1=var(Ydep-MatAGE1*OLSnew_e0);

    
Vect1=Ydep-MatAGE1*Resqinit_e0(:,1);
Vect2=Ydep-MatAGE1*Resqinit_e0(:,Ntau);
b1_e0=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
bL_e0=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));

% Initial conditions: epsilon given AGE
Resqinit_eps=zeros((K3+1)*(K3t+1),Ntau);
Mat=[];
for kk=0:2
   Mat=[Mat hermite(kk,AGE_tot)]; 
end
b_ols = pinv(Mat)*Y_tot;

Ydep_temp=Y_tot - Mat*b_ols + randn(size(Y_tot,1),1);
for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(MatAGE_tot,Ydep_temp,tau);
    Resqinit_eps(:,jtau)=beta1;
end
OLSnew_eps=pinv(MatAGE_tot)*Ydep_temp;
veps=var(Ydep_temp-MatAGE_tot*OLSnew_eps);


Vect1=Ydep_temp-MatAGE_tot*Resqinit_eps(:,1);
Vect2=Ydep_temp-MatAGE_tot*Resqinit_eps(:,Ntau);
b1_eps=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
bL_eps=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    


% % Optional: perturbation of initial conditions
% Resqinit=Resqinit+.1*0.005*randn((K1+1)*(K2+1),Ntau);
% Resqinit_e0=Resqinit_e0+.1*0.005*randn((K3+1)*(K3t+1),Ntau);
% Resqinit_eps=Resqinit_eps+.1*0.005*randn((K4+1)*(K5+1),Ntau);

end