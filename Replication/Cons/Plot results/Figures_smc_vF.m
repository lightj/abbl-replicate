% ABBL (2023)
% Code to generate the figures for the baseline model without heterogeneity

cd '/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Quantile model/'

clear all
clc;
close all;

global Vect Vect_dep xx bdw tau T N K1 K2 K3 K4 Ntau Vectau tau ...
    D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
    Matdraw1 AGE1 Y1 MatAGE1 Matdraw_t Matdraw_lag...
    AGE_t Matdraw_t_lag Y_tot Matdraw_tot MatAGE_tot Y_t Y_lag...
    AGE_tot

% variable for dynamic saving of files
model = 1;

%%
% 1. Load results and extra data for plotting
%%%%%%%%%%%%

% load('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Quantile model/Results/REVISION_cons_smc_l.mat');

% Load the consumption data
Resqfinal_cons=zeros(size(Resqinit_cons));
for jtau=1:Ntau
    for p=1:size(Resqfinal_cons,1)
        Resqfinal_cons(p,jtau)=mean(Resqnew_cons(p,jtau,(2*maxiter/4):maxiter));
    end
end

Resqfinal_a=zeros(size(Resqinit_a));
for jtau=1:Ntau
    for p=1:size(Resqfinal_a,1)
        Resqfinal_a(p,jtau)=mean(Resqnew_a(p,jtau,(2*maxiter/4):maxiter));
    end
end

Resqfinal_a1=zeros(size(Resqinit_a1));
for jtau=1:Ntau
    for p=1:size(Resqfinal_a1,1)
        Resqfinal_a1(p,jtau)=mean(Resqnew_a1(p,jtau,(2*maxiter/4):maxiter));
    end
end

b1_c=mean(mat_b((2*maxiter/4):maxiter,1))
bL_c=mean(mat_b((2*maxiter/4):maxiter,2))
b1_a=mean(mat_b((2*maxiter/4):maxiter,5))
bL_a=mean(mat_b((2*maxiter/4):maxiter,6))
b1_a1=mean(mat_b((2*maxiter/4):maxiter,7))
bL_a1=mean(mat_b((2*maxiter/4):maxiter,8))

Resqtrue=Resqfinal;
Resqtrue_e0=Resqfinal_e0;
Resqtrue_eps=Resqfinal_eps;
Resqtrue_c=Resqfinal_cons;
Resqtrue_a=Resqfinal_a;
Resqtrue_a1=Resqfinal_a1;

b1true=b1;
bLtrue=bL;
b1true_e0=b1_e0;
bLtrue_e0=bL_e0;
b1true_eps=b1_eps;
bLtrue_eps=bL_eps;

b1true_c=b1_c;
bLtrue_c=bL_c;
b1true_a=b1_a;
bLtrue_a=bL_a;
b1true_a1=b1_a1;
bLtrue_a1=bL_a1;

% Load the consumption data
Resqinit= Resqfinal;
Resqinit_e0 = Resqfinal_e0;
Resqinit_eps=Resqfinal_eps;
Resqinit_cons=Resqfinal_cons;
Resqinit_a=Resqfinal_a;
Resqinit_a1=Resqfinal_a1;

%%
% 2. Posterior draws
%%%%%%%%%%%%

Matdraw_final=zeros(N,T);
lik_iter=zeros(N,1);
ESS_iter=zeros(N,T);

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
            
            i_Vect_Matdraw1=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),ua1(:,1),'Uniform',0);
            i_Vect_AGE1=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua1(:,2),'Uniform',0);
            i_Vect_YB=arrayfun(@(x) hermite(x,(YB(iii)-meanYB)/stdYB),ua1(:,3),'Uniform',0);
            i_Vect_ED=arrayfun(@(x) hermite(x,(EDUC(iii)-meanEDUC)/stdEDUC),ua1(:,4),'Uniform',0);
            i_MatA =cat(2,i_Vect_Matdraw1{:}).*cat(2,i_Vect_AGE1{:}).*cat(2,i_Vect_YB{:}).*cat(2,i_Vect_ED{:});
            
            
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
            i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
            i_Vect_Eps_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-particle_draw(:,ttt)-meanY)/stdY),uc(:,3),'Uniform',0);
            i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:});
            
            i_Vect_AGE_t=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
            i_Vect_A_lag=arrayfun(@(x) hermite(x,(i_A(:,ttt-1)-meanA)/stdA),ua(:,1),'Uniform',0);
            i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(particle_draw(:,ttt-1)-meanY)/stdY),ua(:,2),'Uniform',0);
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

Matdraw=Matdraw_final;

Mateta_1=zeros(N,1);
for iii=1:N
  Mateta_1(iii)=Matdraw(iii,ENT(iii));  
end
Mata_true=A;
Matc_true=C;
Mateta_true=Matdraw(:,1:T);
Mateta_1=Matdraw1;
Mateps_true=Y-Mateta_true;
Ytilde=Y;

%%
% 3. Main figures
%%%%%%%%%%%%

%%% FIGURE 4c - Average Insurance

sample=find(D(:));
Vec1=quantile(Mata_true(sample),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

Mat_insurance=zeros(Ntau,Ntau);

sample=find(D);
Mateta_tot =Mateta_true(sample);
Mateps_tot=Mateps_true(sample);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        Vect_AGE_tot=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_tot=arrayfun(@(x) hermite(x,(Vec1(jtau1)-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_tot=arrayfun(@(x) x*hermite(x-1,(Mateta_tot-meanY)/stdY)./stdY,uc(:,2),'Uniform',0);
        Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_tot-meanY)/stdY),uc(:,3),'Uniform',0);
        XX1 = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:});        
        Mat_insurance(jtau1,jtau2)=mean(mean(XX1*Resqtrue_c)); 
    end
end

figs(4)=figure;
set(figs(4), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1))
text(0.8,0.01, 0.85, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 0.75, ['\sigma = ' num2str(std(Mat_insurance(:)),1)])
caxis([0 1])
view([140 15]);
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_c_single_%d.eps',model),'-depsc');
