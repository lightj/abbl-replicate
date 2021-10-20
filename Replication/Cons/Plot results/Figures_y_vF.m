% ABBL (2021)
% Code to generate the figures for the model with heterogeneity but no
% filtering

cd '/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/Consumption/'

clear all
clc;
close all;

% close all

global Vect Vect_dep xx bdw tau T N K1 K2 K3 K4 Ntau Vectau tau ...
    D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
    Matdraw1 AGE1 Y1 MatAGE1 Matdraw_t Matdraw_lag...
    AGE_t Matdraw_t_lag Y_tot Matdraw_tot MatAGE_tot Y_t Y_lag...
    AGE_tot EDUC YB K5 K6 Ybar meanYbar stdYbar

% variable for dynamic saving of files
model = 1;

%%
% 1. Load results and extra data for plotting
%%%%%%%%%%%

load('211015_pmcmc_y.mat');

figs(1)=figure;
plot(mat_lik)
xlim([0, maxiter])
ylabel('Likelihood','FontSize',9)
xlabel('Iteration','FontSize',9)
print(figs(1),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/Consumption/abb_lik_%d.eps',model),'-depsc');

Resqtrue=Resqfinal;
Resqtrue_e0=Resqfinal_e0;
Resqtrue_eps=Resqfinal_eps;
Resqtrue_xi=Resqfinal_xi;
Resqtrue_c=Resqfinal_cons;
Resqtrue_a=Resqfinal_a;
Resqtrue_a1=Resqfinal_a1;

Resqinit_eps=Resqfinal_eps;
Resqinit_xi=Resqfinal_xi;
Resqinit_c=Resqfinal_cons;
Resqinit_a=Resqfinal_a;
Resqinit_a1=Resqfinal_a1;

b1true=b1;
bLtrue=bL;
b1true_e0=b1_e0;
bLtrue_e0=bL_e0;
b1true_eps=b1_eps;
bLtrue_eps=bL_eps;
b1true_xi=b1_xi;
bLtrue_xi=bL_xi;
b1true_c=b1_c;
bLtrue_c=bL_c;
b1true_a=b1_a;
bLtrue_a=bL_a;
b1true_a1=b1_a1;
bLtrue_a1=bL_a1;

% seeds
rng('shuffle')

Mateta_true=zeros(N,T);
Mateps_true=zeros(N,T);
Matc_true=zeros(N,T);
Mata_true=zeros(N,T);
Mateta_1=zeros(N,1);
Matxi_true = zeros(N,1);
Ytilde=zeros(N,T);

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
            end
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
Mata_true=A;
Matc_true=C;
Mateps_true=Y-Mateta_true;
Ytilde=Y;
Matxi_true=Matdraw(:,T+1);

%%
% 3. Main results
%%%%%%%%%%%%

%%% FIGURE 5b - Average Insurance

sample=find(D(:));
Vec1=quantile(Mata_true(sample),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

Mat_insurance=zeros(Ntau,Ntau);

Matxi_tot= Matxi_true;
for tt=2:T
    Matxi_tot=[Matxi_tot; Matxi_true];
end
sample=find(D);
Mateta_tot =Mateta_true(sample);
Mateps_tot=Mateps_true(sample);
Matxi_tot=Matxi_tot(sample);
for jtau1=1:Ntau
    for jtau2=1:Ntau
        Vect_AGE_tot=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_tot=arrayfun(@(x) hermite(x,(Vec1(jtau1)-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_tot=arrayfun(@(x) x*hermite(x-1,(Mateta_tot-meanY)/stdY)./stdY,uc(:,2),'Uniform',0);
        Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_tot-meanY)/stdY),uc(:,3),'Uniform',0);
        Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matxi_tot-meanC)/stdC),uc(:,5),'Uniform',0);
        XX1 = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});        
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
