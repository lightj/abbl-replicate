% ABBL (2021)
% Code to estimate the baseline non-linear consumption model
% with latent earnings but no heterogeneity

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
        Knots_a_eps  Knots_a_a Knots_xi_a Knots_xi_age Knots_xi_eta Knots_a1_eta Knots_a1_age...
        AGE_tot M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13 M14 M15 M16 MatXi MatClag Matdraw MatA MatAlag...
        uc Vect_AGE_tot Vect_A_tot Vect_AGE_t Vect_A_lag Vect_Matdraw_lag Vect_Y_lag Vect_C_lag ua ua1 Vect_AGE1...
        Vect_Matdraw1 meanexpY stdexpY meanexpC stdexpC unit_assets min_val max_val A_tot_RHS;    
    
M1=2;
M2=2;
M3=1;
M4=1;
M5=1;

M6=2;
M7=1;
M8=0;
M9=1;
M10=1;
M11=0;

M12=2;

M13=1;
M14=1;
M15=1;


tol_dnom=2; % Resampling tolerance
maxiter=50; % Stochastic EM iterations
draws=40; % Length of chains in MH step
sigma=1; % S.d. used in initial conditions
Ktime=2; % Only needed if time effects used ( 1 ==> 2003, 2 ==> recession)
MH_scale=5; % Random walk proposal variance is MH_scale*sigma_xi
asset_rule=1;

% Load earnings parameters
load ('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/211015_smc_n.mat')

Nis_scale=5000; % Scaling parameter for importance samples

maxiter=100;
% Load the consumption data
load_cons_data(data);
[Resqinit_cons,Resqinit_a,Resqinit_a1,OLSnew_cons,OLSnew_a,OLSnew_a1,vcons,va,va1,b1_c,bL_c,b1_a,bL_a,b1_a1,bL_a1]=init_cons_smc_abb();

Resqnew_cons=zeros((M1+1)*(M2+1)*(M3+1)*(M4+1),Ntau,maxiter);
Resqnew_a=zeros((M6+1)*(M7+1)*(M8+1)*(M9+1)*(M10+1),Ntau,maxiter);
Resqnew_a1=zeros((M13+1)*(M14+1),Ntau,maxiter);

mat_b=zeros(maxiter,8);

Matdraw=zeros(N,T);
vstore=[];

for iter=1:maxiter
    
    iter
    
    %%%%%%%%%%%%%%% PMCMC %%%%%%%%%%%%%%%
    
    Matdraw_final=zeros(N,T);
    acceptrate=zeros(N,draws);
    lik_iter=zeros(N,1);
    ESS_iter=zeros(N,T);
    tic
    parfor iii=1:N
       
        
        ESS_i=zeros(1,T);
        
        T1=ENT(iii);
        TT=EXT(iii);
        
        Nis=Nis_scale*(TT-T1+1); % No. of importance samples in E-step
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
            i_MatXlag=[i_MatXlag hermite(kk12,(AGE1 -meanAGE)/stdAGE) ];
        end
        
          
                       
            
            i_MatAGE_tot=[];
            for kk3=0:K3
                for kk3t=0:K3t
                    i_MatAGE_tot=[i_MatAGE_tot fun_hermite(kk3,(i_AGE(:)-meanAGE)/stdAGE)...
                        .*fun_hermite(kk3t,i_time(:))];
                end
            end
            
                      
           particle_draw=zeros(Nis,T);

            
            ttt=T1
            while ttt<= TT
                if ttt==T1
                    
                                      
                    
                    
                    alpha=veta1/(veta1 + veps);
                    Mvar=2/((1/veta1)+(1/veps));
                    
                    % Mmu=alpha*(i_Y(1,ttt));
                    
                    Mmu=(1-alpha)*(i_MatAGE1*OLSnew_e0)...
                        + alpha*(i_Y(1,ttt));
                    
                    % Generate the importance samples and their weights
                    particle_draw(:,ttt) = Mmu + mvnrnd(0,Mvar,Nis);
                    
                    % Genreate regressors
                    
                    
                    i_Vect_AGE_tot=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
                    
                    i_Vect_A_tot=arrayfun(@(x) hermite(x,(i_A(:,ttt)-meanA)/stdA),uc(:,1),'Uniform',0);
                    i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
                    i_Vect_Eps_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-particle_draw(:,ttt)-meanY)/stdY),uc(:,3),'Uniform',0);
                    i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:});
                    
                    i_Vect_AGE1=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua1(:,2),'Uniform',0);
                    i_Vect_Matdraw1=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),ua1(:,1),'Uniform',0);
                    i_MatA =cat(2,i_Vect_Matdraw1{:}).*cat(2,i_Vect_AGE1{:});
                    
                    
                    weights=IS_weight_c_a(i_Y,i_C,particle_draw,i_MatAGE_tot,i_MatAGE1,ttt,...
                        Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,Resqinit_cons,...
                        b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,b1_c,bL_c,...
                        Nis,T1,[],i_MatC,Resqinit_a,b1_a,bL_a,Resqinit_a1,b1_a1,bL_a1,i_MatA,i_A,asset_rule);
                    
                        
                    
                    weights=weights./mvnpdf(particle_draw(:,ttt)-Mmu,0,Mvar);
                    %                 weights(isnan(weights))=0;
                    lik=mean(weights);
                    
                elseif ttt>T1
                    
                    alpha=veta/(veta + veps);
                    Mvar=2/((1/veta)+(1/veps));
                    
                    i_Matdraw_lag=[];
                    for kk1=0:K1
                        for kk2=0:K2
                            i_Matdraw_lag=[i_Matdraw_lag hermite(kk1,(particle_draw(:,ttt-1)-meanY)/stdY)...
                                .*hermite(kk2,(i_AGE(:,ttt)-meanAGE)/stdAGE)];
                        end
                    end
                    
                    Mmu=(1-alpha)*(i_Matdraw_lag*OLSnew)...
                        + alpha.*(i_Y(:,ttt));
                                                           
                    % Generate the importance samples and their weights
                    particle_draw(:,ttt) = Mmu + mvnrnd(0,Mvar,Nis);
                    
                    i_Vect_AGE_tot=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
                    i_Vect_A_tot=arrayfun(@(x) hermite(x,(i_A(:,ttt)-meanA)/stdA),uc(:,1),'Uniform',0);
                    i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
                    i_Vect_Eps_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-particle_draw(:,ttt)-meanY)/stdY),uc(:,3),'Uniform',0);
                    i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:});
                    
                    i_Vect_AGE_t=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
                    i_Vect_A_lag=arrayfun(@(x) hermite(x,(i_A(:,ttt-1)-meanA)/stdA),ua(:,1),'Uniform',0);
                    i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(i_Y(:,ttt-1)-meanY)/stdY),ua(:,2),'Uniform',0);
                    i_Vect_Y_lag=arrayfun(@(x) hermite(x,(i_Y(:,ttt-1)-meanY)/stdY),ua(:,3),'Uniform',0);
                    i_Vect_C_lag=arrayfun(@(x) hermite(x,(i_C(:,ttt-1)-meanC)/stdC),ua(:,5),'Uniform',0);
                    i_MatA = cat(2,i_Vect_A_lag{:}).*cat(2,i_Vect_C_lag{:}).*cat(2,i_Vect_Y_lag{:}).*cat(2,i_Vect_Matdraw_lag{:}).*cat(2,i_Vect_AGE_t{:});
                    
                    weights=weights...
                        .*IS_weight_c_a(i_Y,i_C,particle_draw,i_MatAGE_tot,i_MatAGE1,ttt,...
                        Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,Resqinit_cons,...
                        b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,b1_c,bL_c,...
                        Nis,T1,i_Matdraw_lag,i_MatC,Resqinit_a,b1_a,bL_a,Resqinit_a1,b1_a1,bL_a1,i_MatA,i_A,asset_rule);
                    weights=weights./mvnpdf(particle_draw(:,ttt)-Mmu,0,Mvar);
                    
                    % weights(isnan(weights))=0;
                    lik=lik*mean(weights);
                    
                end
                
                % Self-normalization
                weights=weights./sum(weights);
                
                % Effective sample size
                ESS= 1/sum(weights.^2);
                ESS_i(ttt)=ESS;
                %
                % Adaptive resampling (tolerance set by tol)
                if ttt<TT && ESS<tol && ~isnan(ESS)
                    particle_draw = datasample(particle_draw,Nis,'weights',weights);
                    weights=ones(Nis,1)*(1/Nis);
                elseif ttt==TT && ~isnan(ESS)
                    Matdraw_new=datasample(particle_draw,1,'weights',weights);
                elseif isnan(ESS)
                   Matdraw_new=particle_draw(1,:);
                    lik=0.01;
                    ttt=TT;
                end
                ttt=ttt+1;
            end
            
            
              Matdraw_old = Matdraw_new;
            
        
        
        Matdraw_big = Matdraw_old;
        Matdraw_final(iii,:) = Matdraw_big;
        lik_iter(iii)=lik;
%         
         ESS_iter(iii,:)=ESS_i;
%     
    end
    
    toc
    ESS_disp=zeros(1,T);
    ESS_disp2=zeros(1,T);
    for tt=1:T
        ESS_temp=ESS_iter(:,tt);
       sample=find(ESS_temp);
       ESS_disp(tt)=mean( ESS_temp(sample));
       ESS_disp2(tt)=nanmean(ESS_temp(sample));
    end
    ESS_disp    
    ESS_disp2
    
     ESS_store(:,:,iter)=ESS_iter;

    Matdraw=Matdraw_final;
    
        
    options.Display ='off';
    warning off
        
    %%% Generate regressors
     sample=find(D(:));
    
    Matdraw_tot=Matdraw(:,1:T);
    Matdraw_tot=Matdraw_tot(:);
    Matdraw_tot=Matdraw_tot(sample);
    
    Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(Matdraw_tot-meanY)/stdY),uc(:,2),'Uniform',0);
    Vect_Eps_tot=arrayfun(@(x) hermite(x,(Y_tot-Matdraw_tot-meanY)/stdY),uc(:,3),'Uniform',0);
    
    MatClag = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:});
    
    
    
    Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(Y_lag-meanY)/stdY),ua(:,2),'Uniform',0);
    MatAlag = cat(2,Vect_A_lag{:}).*cat(2,Vect_C_lag{:}).*cat(2,Vect_Y_lag{:}).*cat(2,Vect_Matdraw_lag{:}).*cat(2,Vect_AGE_t{:});
    
    
    Matdraw1=zeros(N,1);
    for iii=1:N
    Matdraw1(iii)=Matdraw(iii,ENT(iii));    
    end
    
    Vect_Matdraw1=arrayfun(@(x) hermite(x,(Matdraw1-meanY)/stdY),ua1(:,1),'Uniform',0);
    
    MatA =cat(2,Vect_Matdraw1{:}).*cat(2,Vect_AGE1{:});
    

    
   for jtau=1:Ntau   
        tau=Vectau(jtau);
        Resqnew_cons(:,jtau,iter)=fminunc(@wqregk_cons_age,Resqinit_cons(:,jtau),options);
        if asset_rule==1
        Resqnew_a(:,jtau,iter)=fminunc(@wqregk_a_age,Resqinit_a(:,jtau),options);
        Resqnew_a1(:,jtau,iter)=fminunc(@wqregk_a1_age,Resqinit_a1(:,jtau),options);
        end
    end
    
   
    % Taste heterogeneity: Laplace parameters
    
    Vect1=C_tot-MatClag*Resqnew_cons(:,1,iter);
    Vect2=C_tot-MatClag*Resqnew_cons(:,Ntau,iter);
    b1_c=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_c=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    Resqinit_cons=Resqnew_cons(:,:,iter)

    if asset_rule==1
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
    end
    
    mat_lik(iter)=mean(log(lik_iter));
    mat_lik(iter)
    
    mat_b(iter,1)=b1_c;
    mat_b(iter,2)=bL_c;
%     mat_b(iter,3)=b1_xi;
%     mat_b(iter,4)=bL_xi;
    mat_b(iter,5)=b1_a;
    mat_b(iter,6)=bL_a;
    mat_b(iter,7)=b1_a1;
    mat_b(iter,8)=bL_a1;
    Matdraw;
end

Resqfinal_cons=zeros(size(Resqinit_cons));
for jtau=1:Ntau
    for p=1:size(Resqfinal_cons,1)
        Resqfinal_cons(p,jtau)=mean(Resqnew_cons(p,jtau,(maxiter/2):maxiter));
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

b1_a=mean(mat_b((maxiter/2):maxiter,5))
bL_a=mean(mat_b((maxiter/2):maxiter,6))

b1_a1=mean(mat_b((maxiter/2):maxiter,7))
bL_a1=mean(mat_b((maxiter/2):maxiter,8))

% OLSnew_xi=mean(OLSstore_xi(:,(maxiter/2):maxiter),2);

save 211015_smc_n.mat
