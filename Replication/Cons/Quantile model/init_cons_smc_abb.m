function [Resqinit_cons,Resqinit_a,Resqinit_a1,OLSnew_cons,OLSnew_a,OLSnew_a1,vcons,va,va1,b1_c,bL_c,b1_a,bL_a,b1_a1,bL_a1]=init_cons_smc_abb()

global N K1 K2 K3 K3t K4 K5 Ntau Vectau tau ...
       Y AGE meanAGE stdAGE meanY stdY ENT EXT ...
       MatAGE1 AGE_t MatAGE_tot Y_t Y_lag Y_tot...
       Y1 AGE_tot M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13 M14 M15 M16 M17...
       Matdraw_tot A_tot C_tot A_lag C_lag meanA stdA...
       meanC stdC Matdraw_lag A1 Matdraw1 AGE1  C D T D_t A_t MatXi M12b uc Vect_AGE_tot Vect_A_tot...
       Vect_AGE_t Vect_A_lag Vect_Matdraw_lag Vect_Y_lag Vect_C_lag ua ua1 Vect_AGE1 Vect_Matdraw1...
       meanexpY stdexpY meanexpC stdexpC unit_assets min_val max_val A_tot_RHS Vect_YB Vect_ED...
       YB EDUC meanEDUC stdEDUC meanYB stdYB
   

   
 i_mean_C = zeros(N,1);
 for iii=1:N
     i_mean_C(iii)=mean(C(iii,ENT(iii):EXT(iii)));
 end
 
% Initial conditions: consumption on assets, eta, epsilon, age and xi
sample = find(D);
 i_mean_C_tot = i_mean_C;
 for tt=2:T
     i_mean_C_tot =[i_mean_C_tot;i_mean_C];
 end
 i_mean_C_tot=i_mean_C_tot(sample);
 
Resqinit_cons=zeros((M1+1)*(M2+1)*(M3+1)*(M4+1),Ntau);

[uc1,vc,wc,xc] = ndgrid(0:3);
uc=[uc1(:),vc(:),wc(:),xc(:)];
uc= uc(uc(:,1)<=M1,:);
uc= uc(uc(:,2)<=M2,:);
uc= uc(uc(:,3)<=M3,:);
uc= uc(uc(:,4)<=M4,:);


Vect_A_tot=arrayfun(@(x) hermite(x,(A_tot-meanA)/stdA),uc(:,1),'Uniform',0);
Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(Matdraw_tot-meanY)/stdY),uc(:,2),'Uniform',0);
Vect_Eps_tot=arrayfun(@(x) hermite(x,(Y_tot-Matdraw_tot-meanY)/stdY),uc(:,3),'Uniform',0);
Vect_AGE_tot=arrayfun(@(x) hermite(x,(AGE_tot-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
MatClag = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:});


Cdep=C_tot +randn(size(C_tot,1),1);
for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(MatClag,Cdep,tau);
    Resqinit_cons(:,jtau)=beta1;
end
OLSnew_cons=pinv(MatClag)*Cdep;
vcons=var(Cdep-MatClag*OLSnew_cons);

Vect1=C_tot-MatClag*Resqinit_cons(:,1);
Vect2=C_tot-MatClag*Resqinit_cons(:,Ntau);
b1_c=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
bL_c=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));

% Initial conditions: assets on lag assets, lag consumption, lag Y, lag eta, xi and age

sample = find(D_t);
 i_mean_C_t = i_mean_C;
 for tt=3:T
     i_mean_C_t =[i_mean_C_t;i_mean_C];
 end
 i_mean_C_t=i_mean_C_t(sample);
 
Resqinit_a=zeros((M6+1)*(M7+1)*(M8+1)*(M9+1)*(M10+1),Ntau);

[uc1,vc,wc,xc,yc] = ndgrid(0:3);
ua=[uc1(:),vc(:),wc(:),xc(:),yc(:)];
ua= ua(ua(:,1)<=M6,:);
ua= ua(ua(:,2)<=M7,:);
ua= ua(ua(:,3)<=M8,:);
ua= ua(ua(:,4)<=M9,:);
ua= ua(ua(:,5)<=M10,:);

meanexpY=mean(exp(Y_tot));
stdexpY=std(exp(Y_tot));
meanexpC=mean(exp(C_tot));
stdexpC=std(exp(C_tot));


Vect_A_lag=arrayfun(@(x) hermite(x,(A_lag-meanA)/stdA),ua(:,1),'Uniform',0);
Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(Y_lag - meanY)/stdY),ua(:,2),'Uniform',0);
Vect_Y_lag=arrayfun(@(x) hermite(x,(Y_lag-meanY)/stdY),ua(:,3),'Uniform',0);
Vect_AGE_t=arrayfun(@(x) hermite(x,(AGE_t-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
Vect_C_lag=arrayfun(@(x) hermite(x,(C_lag-meanC)/stdC),ua(:,5),'Uniform',0);

MatAlag = cat(2,Vect_A_lag{:}).*cat(2,Vect_C_lag{:}).*cat(2,Vect_Y_lag{:}).*cat(2,Vect_Matdraw_lag{:}).*cat(2,Vect_AGE_t{:});


Adep=A_t +randn(size(A_t,1),1);
for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(MatAlag,Adep,tau);
    Resqinit_a(:,jtau)=beta1;
end
OLSnew_a=pinv(MatAlag)*Adep;
va=var(Adep-MatAlag*OLSnew_a);

Vect1=A_t-MatAlag*Resqinit_a(:,1);
Vect2=A_t-MatAlag*Resqinit_a(:,Ntau);
b1_a=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
bL_a=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));

% Initial conditions: initial assets
Resqinit_a1=zeros((M13+1)*(M14+1)*(M16+1)*(M17+1),Ntau);

[uc1,vc,xc,yc] = ndgrid(0:1);
ua1=[uc1(:),vc(:),xc(:),yc(:)];


Vect_Matdraw1=arrayfun(@(x) hermite(x,(Matdraw1-meanY)/stdY),ua1(:,1),'Uniform',0);
Vect_AGE1=arrayfun(@(x) hermite(x,(AGE1-meanAGE)/stdAGE),ua1(:,2),'Uniform',0);
Vect_YB=arrayfun(@(x) hermite(x,(YB-meanYB)/stdYB),ua1(:,3),'Uniform',0);
Vect_ED=arrayfun(@(x) hermite(x,(EDUC-meanEDUC)/stdEDUC),ua1(:,4),'Uniform',0);

MatA =cat(2,Vect_Matdraw1{:}).*cat(2,Vect_AGE1{:}).*cat(2,Vect_YB{:}).*cat(2,Vect_ED{:});


for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(MatA,A1,tau);
    Resqinit_a1(:,jtau)=beta1;
end
OLSnew_a1=pinv(MatA)*(A1);
va1=var(A1-MatA*OLSnew_a1);

Vect1=A1-MatA*Resqinit_a1(:,1);
Vect2=A1-MatA*Resqinit_a1(:,Ntau);
b1_a1=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
bL_a1=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));


end