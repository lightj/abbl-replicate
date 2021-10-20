% ABBL (2021)
% Code for model with heterogeneity and not filtering

clear all
clc;
cd '/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/Consumption/'

global T N K1 K2 K3 K4 K5 Ntau Vectau tau AGE1 Y1 A1...
    D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
    Matdraw1 MatAGE1 Matdraw_t Matdraw_lag AGE_t Mat_ETA_AGE Y_tot Matdraw_tot MatAGE_tot...
    Y_t Y_lag meanT stdT K3t Matdraw_eps C A A_tot  C_tot meanA stdA meanC stdC A_t A_lag C_lag ...
    Vect_c_xi  Vect_c_age Vect_c_eta Vect_c_eps Vect_c_a  Vect_a_xi Vect_a_age Vect_a_eta ...
    Vect_a_eps Vect_a_c Vect_a_a Vect_xi_a Vect_xi_age Vect_xi_eta Vect_a1_eta Vect_a1_age...
    Knots_c_xi Knots_c_age Knots_c_eta Knots_c_eps Knots_c_a Knots_a_xi Knots_a_age Knots_a_eta ...
    Knots_a_eps Knots_a_c Knots_a_a Knots_xi_a Knots_xi_age Knots_xi_eta Knots_a1_eta Knots_a1_age...
    AGE_tot M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13 M14 M15 M16 MatXi MatClag Matdraw MatA MatAlag...
    uc Vect_AGE_tot Vect_A_tot Vect_AGE_t Vect_A_lag Vect_Matdraw_lag Vect_Y_lag Vect_C_lag ua ua1 Vect_AGE1 Vect_Matdraw1 unit_interval...
    EDUC YB meanEDUC meanYB stdEDUC stdYB Ybar  meanYbar stdYbar M12a M12b M17   Vect_YB Vect_ED...
    Vect_Matdraw_tot

M1=2;
M2=2;
M3=0;
M4=1;
M5=1;
K=2;
M6=2;
M7=1;
M8=0;
M9=1;
M10=1;
M11=1;

M12=1;
M12a=1;
M12b=1;

M13=1;
M14=1;
M15=1;
M16=1;
M17=1;

load ('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/211015_smc_n.mat')

asset_rule=1;
Nis_scale=1; % Scaling parameter for importance samples
tol_dnom=2; % Resampling tolerance
draws=50; % Length of chains in MH step
sigma=1; % S.d. used in initial conditions
Ktime=2; % Only needed if time effects used ( 1 ==> 2003, 2 ==> recession)
MH_scale=5; % Random walk proposal variance is MH_scale*sigma_xi
unit_interval=0;
maxiter=1600;

% Load the consumption data
load_cons_data(data);
[Resqinit_cons,Resqinit_a,Resqinit_xi,Resqinit_a1,OLSnew_cons,OLSnew_a,OLSnew_xi,OLSnew_a1,vcons,va,vxi,va1,b1_c,bL_c,b1_a,bL_a,b1_a1,bL_a1,b1_xi,bL_xi]=init_cons_pmcmc_y();
Resqinit= Resqfinal;
Resqinit_e0 = Resqfinal_e0;
Resqinit_eps=Resqfinal_eps;

Resqnew_cons=zeros((M1+1)*(M2+1)*(M3+1)*(M4+1)*(M5+1),Ntau,maxiter);
Resqnew_a=zeros((M6+1)*(M7+1)*(M8+1)*(M9+1)*(M10+1)*(M11+1),Ntau,maxiter);
Resqnew_xi=zeros((M12+1)*(M12a+1)*(M12b+1),Ntau,maxiter);
Resqnew_a1=zeros((M13+1)*(M14+1)*(M15+1)*(M16+1)*(M17+1),Ntau,maxiter);

% b1_c=5;
% bL_c=5;
% b1_a=5;
% bL_a=5;
% b1_xi=5;
% bL_xi=5;
% b1_a1=5;
% bL_a1=5;

mat_b=zeros(maxiter,8);

Matdraw=zeros(N,T+1);
vstore=[];
xi_draws=zeros(N,1);
ESS_store=zeros(N,T,maxiter);

mat_lik=zeros(maxiter,1);

for iter=1:maxiter
    
    iter
    
    %%%%%%%%%%%%%%% MCMC %%%%%%%%%%%%%%%
    
    Matdraw_final=zeros(N,T+1);
    acceptrate=zeros(N,draws);
    lik_iter=zeros(N,1);
    ESS_iter=zeros(N,T);
    
    parfor iii=1:N
        
        a_i_old=[];
        Matdraw_old=[];
        acc=zeros(1,draws);
        
        ESS_i=zeros(draws,T);
        lik_i=zeros(draws,1);
        
        T1=ENT(iii);
        TT=EXT(iii);
        
        Nis=Nis_scale; % No. of importance samples in E-step
        tol=Nis/tol_dnom;
        
        i_Y=ones(Nis,1)*Y(iii,:);
        i_AGE=ones(Nis,1)*AGE(iii,:);
        
        i_A=ones(Nis,1)*A(iii,:);
        i_C=ones(Nis,1)*C(iii,:);
        i_time=0;
        
        i_MatAGE1=[];
        for kk4=0:K4
            for kk5=0:K5
                for kk6=0:K6
                    i_MatAGE1=[i_MatAGE1 hermite(kk4,(i_AGE(:,T1)-meanAGE)/stdAGE)...
                        .*hermite(kk5,(EDUC(iii)-meanEDUC)/stdEDUC)...
                        .*hermite(kk6,(YB(iii)-meanYB)/stdYB)];
                end
            end
        end
        
        i_MatXlag=[];
        for kk12=0:M12
            for kk13=0:M12a
                for kk14=0:M12b
                    i_MatXlag=[i_MatXlag hermite(kk12,(EDUC(iii)-meanEDUC)/stdEDUC)...
                        .*hermite(kk13,(YB(iii)-meanYB)/stdYB)...
                        .*hermite(kk14,(Ybar(iii) -meanYbar)/stdYbar)];
                end;
            end
        end
        
        
        j=1;
        while j<=draws
            
            % Define useful matrices
            particle_draw=zeros(Nis,T);
            mat_resamp=zeros(T,1);
            
            
            if j==1
                if iter==1
                    a_i_new =i_MatXlag(1,:)*OLSnew_xi + mvnrnd(0,MH_scale*vxi,1);
                elseif iter>1
                    a_i_new = xi_draws(iii)
                end
            elseif j>1
                a_i_new = a_i_old + mvnrnd(0,MH_scale*vxi,1);
            end
            
            
            i_MatAGE_tot=[];
            for kk3=0:K3
                for kk3t=0:K3t
                    i_MatAGE_tot=[i_MatAGE_tot fun_hermite(kk3,(i_AGE(:)-meanAGE)/stdAGE)...
                        .*fun_hermite(kk3t,i_time(:))];
                end
            end
            
            
            
            
            ttt=T1
            while ttt<= TT
                if ttt==T1
                    
                    i_Vect_AGE_tot=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
                    i_Vect_A_tot=arrayfun(@(x) hermite(x,(i_A(:,ttt)-meanA)/stdA),uc(:,1),'Uniform',0);
                    i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
                    i_Vect_Xi_tot=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),uc(:,5),'Uniform',0);
                    i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_AGE_tot{:}).*cat(2,i_Vect_Xi_tot{:});
                    
                    
                    i_Vect_AGE1=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua1(:,2),'Uniform',0);
                    i_Vect_Matdraw1=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-meanY)/stdY),ua1(:,1),'Uniform',0);
                    i_Vect_Xi1=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),ua1(:,3),'Uniform',0);
                    i_Vect_YB=arrayfun(@(x) hermite(x,(YB(iii)-meanYB)/stdYB),ua1(:,4),'Uniform',0);
                    i_Vect_ED=arrayfun(@(x) hermite(x,(EDUC(iii)-meanEDUC)/stdEDUC),ua1(:,5),'Uniform',0);
                    
                    i_MatA =cat(2,i_Vect_Matdraw1{:}).*cat(2,i_Vect_AGE1{:}).*cat(2,i_Vect_Xi1{:}).*cat(2,i_Vect_YB{:}).*cat(2,i_Vect_ED{:});
                    
                    
                    
                    weights=mh_likl(i_Y,i_C,particle_draw,i_MatAGE_tot,i_MatAGE1,ttt,...
                        Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,Resqinit_cons,...
                        b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,b1_c,bL_c,...
                        Nis,T1,[],i_MatC,Resqinit_a,b1_a,bL_a,Resqinit_a1,b1_a1,bL_a1,i_MatA,i_A,asset_rule);
                    
                    
                elseif ttt>T1
                    
                    
                    i_Vect_AGE_tot=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
                    i_Vect_A_tot=arrayfun(@(x) hermite(x,(i_A(:,ttt)-meanA)/stdA),uc(:,1),'Uniform',0);
                    i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
                    i_Vect_Xi_tot=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),uc(:,5),'Uniform',0);
                    i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_AGE_tot{:}).*cat(2,i_Vect_Xi_tot{:});
                    
                    
                    i_Vect_AGE_t=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
                    i_Vect_A_lag=arrayfun(@(x) hermite(x,(i_A(:,ttt-1)-meanA)/stdA),ua(:,1),'Uniform',0);
                    i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(i_Y(:,ttt-1)-meanY)/stdY),ua(:,2),'Uniform',0);
                    i_Vect_C_lag=arrayfun(@(x) hermite(x,(i_C(:,ttt-1)-meanC)/stdC),ua(:,5),'Uniform',0);
                    i_Vect_Xi_lag=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),ua(:,6),'Uniform',0);
                    
                    i_MatA = cat(2,i_Vect_A_lag{:}).*cat(2,i_Vect_C_lag{:}).*cat(2,i_Vect_Matdraw_lag{:}).*cat(2,i_Vect_AGE_t{:}).*cat(2,i_Vect_Xi_lag{:});
                    
                    
                    weights=weights...
                        .*mh_likl(i_Y,i_C,particle_draw,i_MatAGE_tot,i_MatAGE1,ttt,...
                        Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,Resqinit_cons,...
                        b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,b1_c,bL_c,...
                        Nis,T1,[],i_MatC,Resqinit_a,b1_a,bL_a,Resqinit_a1,b1_a1,bL_a1,i_MatA,i_A,asset_rule);
                    
                end
                ttt=ttt+1;
            end
            
            newObj=weights*fun_prior_c(a_i_new,i_MatXlag,Ntau,Vectau,Resqinit_xi,b1_xi,bL_xi);
            
            
            if j==1
                Matdraw_old = particle_draw;
                a_i_old = a_i_new;
                Obj_chain=[newObj zeros(1,draws-1)];
                lik_chain=[weights zeros(1,draws-1)];

            elseif j>1
                r=(min([1 newObj./Obj_chain(:,j-1)]'))';
                prob=rand(1,1);
                Obj_chain(:,j)=(prob<=r).*newObj+(prob>r).*Obj_chain(:,j-1);
                lik_chain(:,j)=(prob<=r).*weights+(prob>r).*lik_chain(:,j-1);

                Matdraw_old=(prob<=r).*particle_draw+(prob>r).*Matdraw_old;
                a_i_old=(prob<=r).*a_i_new+(prob>r).*a_i_old;
                acc(j)=(prob<=r);
            end
            
            likl=mean(lik_chain(draws/2 : end));
            j=j+1;
            
        end
        
        Matdraw_big = [Matdraw_old a_i_old];
        acceptrate(iii,:)=acc;
        Matdraw_final(iii,:) = Matdraw_big;
        lik_iter(iii)=likl;

    end
     
    Matdraw=Matdraw_final;
    Acceptrate_store(:,:,iter)=acceptrate;
    mean(acceptrate)
    
    options.Display ='off';
    warning off
    
    %% Normalize xi
    if iter > 0
       
        C1=zeros(N,1);
        for iii=1:N
            C1(iii)=C(iii,ENT(iii));
        end
        

        
        Vect_AGE_temp=arrayfun(@(x) hermite(x,(AGE1-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_temp=arrayfun(@(x) hermite(x,(A1-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_temp=arrayfun(@(x) hermite(x,(Y1-meanY)/stdY),uc(:,2),'Uniform',0);
        Vect_Xi_temp=arrayfun(@(x) hermite(x,(Matdraw(:,T+1))),uc(:,5),'Uniform',0);
        
        MatClag = cat(2,Vect_A_temp{:}).*cat(2,Vect_Matdraw_temp{:}).*cat(2,Vect_AGE_temp{:}).*cat(2,Vect_Xi_temp{:});
        
        
        b_ols=pinv(MatClag)*C1;
        
        Vect_AGE_temp=arrayfun(@(x) hermite(x,0),uc(:,4),'Uniform',0);
        Vect_A_temp=arrayfun(@(x) hermite(x,0),uc(:,1),'Uniform',0);
        Vect_Matdraw_temp=arrayfun(@(x) hermite(x,0),uc(:,2),'Uniform',0);
        Vect_Xi_temp=arrayfun(@(x) hermite(x,(Matdraw(:,T+1))),uc(:,5),'Uniform',0);
        
        MatClag = cat(2,Vect_A_temp{:}).*cat(2,Vect_Matdraw_temp{:}).*cat(2,Vect_AGE_temp{:}).*cat(2,Vect_Xi_temp{:});
        
        
        Matdraw(:,T+1) = MatClag*b_ols;
        
    end
    
    %%% Generate regressors
    sample=find(D(:));

    Matzeta_tot=Matdraw(:,T+1);
    for j=2:T
        Matzeta_tot=[Matzeta_tot;Matdraw(:,T+1)];
    end
    Matzeta_tot=Matzeta_tot(sample);
    

    Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matzeta_tot-meanC)/stdC),uc(:,5),'Uniform',0);
    MatClag = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});
    
    
    
    sample=find(D_t);
    i_mean_C_t=Matdraw(:,T+1);
    for tt=3:T
        i_mean_C_t=[i_mean_C_t Matdraw(:,T+1)];
    end
    i_mean_C_t=i_mean_C_t(sample);
 
    
    Vect_Xi_lag=arrayfun(@(x) hermite(x,(i_mean_C_t-meanC)/stdC),ua(:,6),'Uniform',0);
    
    MatAlag = cat(2,Vect_A_lag{:}).*cat(2,Vect_C_lag{:}).*cat(2,Vect_Matdraw_lag{:}).*cat(2,Vect_AGE_t{:}).*cat(2,Vect_Xi_lag{:});
    
    
    Vect_Xi1=arrayfun(@(x) hermite(x,(Matdraw(:,T+1)-meanC)/stdC),ua1(:,3),'Uniform',0);
    MatA =cat(2,Vect_Matdraw1{:}).*cat(2,Vect_AGE1{:}).*cat(2,Vect_Xi1{:}).*cat(2,Vect_YB{:}).*cat(2,Vect_ED{:});
    
    
    for jtau=1:Ntau
        tau=Vectau(jtau);
        Resqnew_cons(:,jtau,iter)=fminunc(@wqregk_cons_age,Resqinit_cons(:,jtau),options);
        Resqnew_xi(:,jtau,iter)=fminunc(@wqregk_xi_age,Resqinit_xi(:,jtau),options);
        Resqnew_a(:,jtau,iter)=fminunc(@wqregk_a_age,Resqinit_a(:,jtau),options);
        Resqnew_a1(:,jtau,iter)=fminunc(@wqregk_a1_age,Resqinit_a1(:,jtau),options);
    end
    
    beta=pinv(MatXi)*Matdraw(:,T+1);
    vxi=var(Matdraw(:,T+1)-MatXi*beta)
    
    %     Resqnew_cons(17,:,iter)=Resqnew_cons(17,:,iter)-mean(Resqnew_cons(17,:,iter)')'*ones(1,Ntau) + stdC;
    %     Resqnew_cons(1,:,iter)=Resqnew_cons(1,:,iter)-mean(Resqnew_cons(1,:,iter)')'*ones(1,Ntau) + meanC;
    %     Resqnew_cons(1,:,iter)=Resqnew_cons(1,:,iter)-((1-Vectau(Ntau))/bL_c-Vectau(1)/b1_c)*ones(1,Ntau);
    %
    %
    % Taste heterogeneity: Laplace parameters
    
    Vect1=C_tot-MatClag*Resqnew_cons(:,1,iter);
    Vect2=C_tot-MatClag*Resqnew_cons(:,Ntau,iter);
    b1_c=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_c=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    Resqinit_cons=Resqnew_cons(:,:,iter)
    
    Vect1=Matdraw(:,T+1)-MatXi*Resqnew_xi(:,1,iter);
    Vect2=Matdraw(:,T+1)-MatXi*Resqnew_xi(:,Ntau,iter);
    b1_xi=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_xi=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    Resqinit_xi=Resqnew_xi(:,:,iter)
    
    Vect1=A_t-MatAlag*Resqnew_a(:,1,iter);
    Vect2=A_t-MatAlag*Resqnew_a(:,Ntau,iter);
    b1_a=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_a=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    
    Vect1=A1-MatA*Resqnew_a1(:,1,iter);
    Vect2=A1-MatA*Resqnew_a1(:,Ntau,iter);
    b1_a1=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_a1=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    
    Resqinit_a=Resqnew_a(:,:,iter)
    Resqinit_a1=Resqnew_a1(:,:,iter)
    
    
    mat_lik(iter)=mean(log(lik_iter));
    mat_lik(iter)
    
    mat_b(iter,1)=b1_c;
    mat_b(iter,2)=bL_c;
    mat_b(iter,3)=b1_xi;
    mat_b(iter,4)=bL_xi;
    mat_b(iter,5)=b1_a;
    mat_b(iter,6)=bL_a;
    mat_b(iter,7)=b1_a1;
    mat_b(iter,8)=bL_a1;
    
    vstore=[vstore vxi];
    xi_draws=Matdraw(:,T+1);
end

Resqfinal_cons=zeros(size(Resqinit_cons));
for jtau=1:Ntau
    for p=1:size(Resqfinal_cons,1)
        Resqfinal_cons(p,jtau)=mean(Resqnew_cons(p,jtau,(maxiter/2):maxiter));
    end
end

Resqfinal_xi=zeros(size(Resqinit_xi));
for jtau=1:Ntau
    for p=1:size(Resqfinal_xi,1)
        Resqfinal_xi(p,jtau)=mean(Resqnew_xi(p,jtau,(maxiter/2):maxiter));
    end
end

Resqfinal_a=zeros(size(Resqinit_a));
for jtau=1:Ntau
    for p=1:size(Resqfinal_a,1)
        Resqfinal_a(p,jtau)=mean(Resqnew_a(p,jtau,(maxiter/2):maxiter));
    end
end

Resqfinal_a1=zeros(size(Resqinit_a1));
for jtau=1:Ntau
    for p=1:size(Resqfinal_a1,1)
        Resqfinal_a1(p,jtau)=mean(Resqnew_a1(p,jtau,(maxiter/2):maxiter));
    end
end

b1_c=mean(mat_b((maxiter/2):maxiter,1))
bL_c=mean(mat_b((maxiter/2):maxiter,2))
b1_xi=mean(mat_b((maxiter/2):maxiter,3))
bL_xi=mean(mat_b((maxiter/2):maxiter,4))
b1_a=mean(mat_b((maxiter/2):maxiter,5))
bL_a=mean(mat_b((maxiter/2):maxiter,6))
b1_a1=mean(mat_b((maxiter/2):maxiter,7))
bL_a1=mean(mat_b((maxiter/2):maxiter,8))

% OLSnew_xi=mean(OLSstore_xi(:,(maxiter/2):maxiter),2);
%vxi=mean(var_store(1,(maxiter/2):maxiter),2);




save 211015_pmcmc_y.mat
