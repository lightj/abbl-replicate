% ABBL (2021)
% Code to generate the figures for the model with asset rule

cd '/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/Consumption/'

clear all
clc;
close all;

global Vect Vect_dep xx bdw tau T N K1 K2 K3 K4 Ntau Vectau tau ...
    D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
    Matdraw1 AGE1 Y1 MatAGE1 Matdraw_t Matdraw_lag...
    AGE_t Matdraw_t_lag Y_tot Matdraw_tot MatAGE_tot Y_t Y_lag...
    AGE_tot EDUC YB K5 K6 Ybar meanYbar stdYbar

% variable for dynamic saving of files
model = 2000;

%%
% 1. Load results and extra data for plotting
%%%%%%%%%%%%

load ('/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/Cons/Quantile Model/211015_pmcmc_n_assets.mat')

% Check marginal likelihood
figs(1)=figure;
plot(mat_lik)
xlim([0, maxiter])
ylabel('Likelihood','FontSize',9)
xlabel('Iteration','FontSize',9)
print(figs(1),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_lik_%d.eps',model),'-depsc');

Resqinit= Resqfinal;
Resqinit_e0 = Resqfinal_e0;
Resqinit_eps=Resqfinal_eps;

Resqfinal_cons=zeros(size(Resqinit_cons));
for jtau=1:Ntau
    for p=1:size(Resqfinal_cons,1)
        Resqfinal_cons(p,jtau)=mean(Resqnew_cons(p,jtau,(2*maxiter/4):maxiter));
    end
end

Resqfinal_xi=zeros(size(Resqinit_xi));
for jtau=1:Ntau
    for p=1:size(Resqfinal_xi,1)
        Resqfinal_xi(p,jtau)=mean(Resqnew_xi(p,jtau,(2*maxiter/4):maxiter));
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
b1_xi=mean(mat_b((2*maxiter/4):maxiter,3))
bL_xi=mean(mat_b((2*maxiter/4):maxiter,4))
b1_a=mean(mat_b((2*maxiter/4):maxiter,5))
bL_a=mean(mat_b((2*maxiter/4):maxiter,6))
b1_a1=mean(mat_b((2*maxiter/4):maxiter,7))
bL_a1=mean(mat_b((2*maxiter/4):maxiter,8))

Resqtrue=Resqfinal;
Resqtrue_e0=Resqfinal_e0;
Resqtrue_eps=Resqfinal_eps;
Resqtrue_xi=Resqfinal_xi;
Resqtrue_c=Resqfinal_cons;
Resqtrue_a=Resqfinal_a;
Resqtrue_a1=Resqfinal_a1;

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

Resqinit_cons=Resqfinal_cons;
Resqinit_xi=Resqfinal_xi;
Resqinit_a=Resqfinal_a;
Resqinit_a1=Resqfinal_a1;

%%
% 2. Simulate posterior draws
%%%%%%%%%%%%

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
                i_Vect_Xi_tot=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),uc(:,5),'Uniform',0);
                i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:}).*cat(2,i_Vect_Xi_tot{:});
                
                
                i_Vect_AGE1=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua1(:,2),'Uniform',0);
                i_Vect_Matdraw1=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),ua1(:,1),'Uniform',0);
                i_Vect_Xi1=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),ua1(:,3),'Uniform',0);
                i_Vect_YB=arrayfun(@(x) hermite(x,(YB(iii)-meanYB)/stdYB),ua1(:,4),'Uniform',0);
                i_Vect_ED=arrayfun(@(x) hermite(x,(EDUC(iii)-meanEDUC)/stdEDUC),ua1(:,5),'Uniform',0);
                
                i_MatA =cat(2,i_Vect_Matdraw1{:}).*cat(2,i_Vect_AGE1{:}).*cat(2,i_Vect_Xi1{:}).*cat(2,i_Vect_YB{:}).*cat(2,i_Vect_ED{:});
                
                
                
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
                        for kk2t=0:K2t
                            i_Matdraw_lag=[i_Matdraw_lag hermite(kk1,(particle_draw(:,ttt-1)-meanY)/stdY)...
                                .*hermite(kk2,(i_AGE(:,ttt)-meanAGE)/stdAGE)...
                                .*hermite(kk2t,ttt)];
                        end
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
                i_Vect_Xi_tot=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),uc(:,5),'Uniform',0);
                i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:}).*cat(2,i_Vect_Xi_tot{:});
                
                
                i_Vect_AGE_t=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
                i_Vect_A_lag=arrayfun(@(x) hermite(x,(i_A(:,ttt-1)-meanA)/stdA),ua(:,1),'Uniform',0);
                i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(particle_draw(:,ttt-1)-meanY)/stdY),ua(:,2),'Uniform',0);
                i_Vect_Y_lag=arrayfun(@(x) hermite(x,(i_Y(:,ttt-1)-particle_draw(:,ttt-1) - meanY)/stdY),ua(:,3),'Uniform',0);
                i_Vect_C_lag=arrayfun(@(x) hermite(x,(i_C(:,ttt-1)-meanC)/stdC),ua(:,5),'Uniform',0);
                i_Vect_Xi_lag=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),ua(:,6),'Uniform',0);
                
                i_MatA = cat(2,i_Vect_A_lag{:}).*cat(2,i_Vect_C_lag{:}).*cat(2,i_Vect_Y_lag{:}).*cat(2,i_Vect_Matdraw_lag{:}).*cat(2,i_Vect_AGE_t{:}).*cat(2,i_Vect_Xi_lag{:});
                
                
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
            ESS_i(j,ttt)=ESS;
            %
            % Adaptive resampling (tolerance set by tol)
            if ttt<TT && ESS<tol && ~isnan(ESS)
                particle_draw = datasample(particle_draw,Nis,'weights',weights);
                weights=ones(Nis,1)*(1/Nis);
            elseif ttt==TT && ~isnan(ESS)
                Matdraw_new=datasample(particle_draw,1,'weights',weights);
            elseif isnan(ESS)
                if j>1
                    Matdraw_new=Matdraw_old;
                else
                    Matdraw_new=particle_draw(1,:);
                end
                lik=0.01;
                ttt=TT;
            end
            ttt=ttt+1;
        end
        
        newObj=lik*fun_prior_c(a_i_new,i_MatXlag,Ntau,Vectau,Resqinit_xi,b1_xi,bL_xi);
        
        if j==1
            Matdraw_old = Matdraw_new;
            a_i_old = a_i_new;
            Obj_chain=[newObj zeros(1,draws-1)];
            lik_chain=lik;
        elseif j>1
            r=(min([1 newObj./Obj_chain(:,j-1)]'))';
            prob=rand(1,1);
            Obj_chain(:,j)=(prob<=r).*newObj+(prob>r).*Obj_chain(:,j-1);
            Matdraw_old=(prob<=r).*Matdraw_new+(prob>r).*Matdraw_old;
            a_i_old=(prob<=r).*a_i_new+(prob>r).*a_i_old;
            lik_chain=(prob<=r).*lik+(prob>r).*lik_chain;
            acc(j)=(prob<=r);
        end
        
        
        j=j+1;
        
    end
    
    Matdraw_big = [Matdraw_old a_i_old];
    acceptrate(iii,:)=acc;
    Matdraw_final(iii,:) = Matdraw_big;
    lik_iter(iii)=lik_chain;
    
    ESS_mean_i=mean(ESS_i);
    ESS_iter(iii,:)=ESS_mean_i;
end

Matdraw=Matdraw_final;
Matdraw1=zeros(N,1);
for iii=1:N
    Matdraw1(iii)=Matdraw(iii,ENT(iii));
end

Mata_true=A;
Matc_true=C;
Mateta_true=Matdraw(:,1:T);
Mateta_1=Matdraw1;
Mateps_true=Y-Mateta_true;
Ytilde=Y;
Matxi_true=Matdraw(:,T+1);


%%% FIGURE 5 - Marginal Distributions
figs(5) = figure;
set(figs(5), 'Position', [10 10 900 300]);

% C
Vect=Matc_true(:,1:T);
Vect=Vect(:);
sample=find(D);
Vect=Vect(sample);
[f_eta,xi]=ksdensity(Vect);
subplot(1,4,1)
plot(xi,f_eta,'-','Linewidth',1,'Color','b')
xlim([-5 5])
hold on
Vect=C(:,1:T);
Vect=Vect(:);
sample=find(D);
Vect=Vect(sample);
[f_eta,xi]=ksdensity(Vect);
plot(xi,f_eta,'-','Linewidth',1,'Color','g')
title('Consumption')

% A
Vect=Mata_true(:,2:T);
Vect=Vect(:);
sample=find(D_t);
Vect=Vect(sample);
[f_eta,xi]=ksdensity(Vect);
subplot(1,4,2)
plot(xi,f_eta,'-','Linewidth',1,'Color','b')
hold on
Vect=A(:,2:T);
Vect=Vect(:);
sample=find(D_t);
Vect=Vect(sample);
[f_eta,xi]=ksdensity(Vect);
plot(xi,f_eta,'-','Linewidth',1,'Color','g')
title('Assets')

% A1
Vect=zeros(N,1);
for ii=1:N
   Vect(ii)=Mata_true(ii,ENT(ii)); 
end
[f_eta,xi]=ksdensity(Vect);
subplot(1,4,3)
plot(xi,f_eta,'-','Linewidth',1,'Color','b')
hold on
Vect=zeros(N,1);
for ii=1:N
   Vect(ii)=A(ii,ENT(ii)); 
end
[f_eta,xi]=ksdensity(Vect);
plot(xi,f_eta,'-','Linewidth',1,'Color','g')
title('Initial Assets')

% A1
Vect=Matxi_true(:);
[f_eta,xi]=ksdensity(Vect);
subplot(1,4,4)
plot(xi,f_eta,'-','Linewidth',1,'Color','b')
title('Heterogeneity')

print(figs(5),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/densities2_%d.eps',model),'-depsc');

%%% FIGURE 5

figs(13) = figure;
set(figs(13), 'Position', [10 10 600 300]);

% ETA
sample=find((EDUC(:)==0).*(YB(:)<=1969));
Vect=Matxi_true(sample);
[f_eta,xi]=ksdensity(Vect);
plot(xi,f_eta,'-','Linewidth',1,'Color','b')
hold on
sample=find((EDUC(:)==1).*(YB(:)<=1969));
Vect=Matxi_true(sample);
[f_eta,xi]=ksdensity(Vect);
plot(xi,f_eta,'-','Linewidth',1,'Color','r')
sample=find((EDUC(:)==0).*(YB(:)>=1969));
Vect=Matxi_true(sample);
[f_eta,xi]=ksdensity(Vect);
plot(xi,f_eta,'-','Linewidth',1,'Color','g')
sample=find((EDUC(:)==1).*(YB(:)>=1969));
Vect=Matxi_true(sample);
[f_eta,xi]=ksdensity(Vect);
plot(xi,f_eta,'-','Linewidth',1,'Color','c')
hold off
xlabel('\xi','FontSize',9)
ylabel('Density','FontSize',9)
title('\xi density estimate','FontSize',9)
legend('Old, No College','Old, College','Young, No College','Young, College')
print(figs(13),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/densities_xi_%d.eps',model),'-depsc');

%%% FIGURE 6 - Main Insurance \eta

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

% Figure 5c
figs(4)=figure;
set(figs(4), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[-0.1 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1))
view([140 15]);
caxis([0 1])
text(0.8,0.01, 0.85, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 0.75, ['\sigma = ' num2str(std(Mat_insurance(:)),2)])
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_c_single_%d.eps',model),'-depsc');

%%% FIGURE 6 - Main Insurance Y
XX=[];
for kk1=0:M1
    for kk2=0:M2
        for kk3=0:M4
            XX=[XX hermite(kk1,(A_tot-meanA)/stdA)...
                        .*hermite(kk2,(Y_tot-meanY)/stdY)...
                        .*hermite(kk3,(AGE_tot-meanAGE)/stdAGE)];
        end
    end
end

Resqdata_cons=zeros((M1+1)*(M2+1)*(M4+1),Ntau);
for jtau=1:Ntau
    beta1=rq(XX,C_tot,Vectau(jtau));
    Resqdata_cons(:,jtau)=beta1;
end

sample=find(D(:));
Vec1=quantile(Mata_true(sample),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

Mat_insurance=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        
        XX=[];
        for kk1=0:M1
            for kk2=0:M2
                for kk3=0:M4
                    if kk2<1
                        XX=[XX zeros(size(Y_tot,1),1)];
                    else
                        XX=[XX hermite(kk1,(Vec1(jtau1)-meanA)/stdA)...
                            .*kk2/stdY*hermite(kk2-1,(Y_tot-meanY)/stdY)...
                            .*hermite(kk3,(Vec2(jtau2)-meanAGE)/stdAGE)];
                    end
                end
            end
        end
        Mat_insurance(jtau1,jtau2)=mean(mean(XX*Resqdata_cons));

    end
end

% Figure 5c
figs(4)=figure;
set(figs(4), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_insurance)
xlabel('percentile \tau_{age}','FontSize',8)
ylabel('percentile \tau_{assets}','FontSize',8)
zlabel('response to \eta','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1))
view([140 15]);
caxis([0 1])
title([' \mu =',num2str(mean(mean(Mat_insurance))),', \sigma =',num2str(std((Mat_insurance(:))))])
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_cy_single_%d.eps',model),'-depsc');

%%%% MAIN HETEROGENEITY FIGURES

Ntau_plot=5;
quant_plot_grid = [0.1 0.25 0.5 0.75 0.9];

sample=find(D(:));
Mateta_tot =Mateta_true(sample);
Mateps_tot=Mateps_true(sample);

Vec1=quantile(Mata_true(sample),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

%%% (a) consumption derivative wrt eta 
for it = 1:5
tau_use = quant_plot_grid(it);
Matxi_use = quantile(Matxi_true,tau_use);

Mat_insurance=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        
        
        Vect_AGE_tot=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_tot=arrayfun(@(x) hermite(x,(Vec1(jtau1)-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_tot=arrayfun(@(x) x*hermite(x-1,(Mateta_tot-meanY)/stdY)./stdY,uc(:,2),'Uniform',0);
        Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_tot-meanY)/stdY),uc(:,3),'Uniform',0);
        Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matxi_use-meanC)/stdC),uc(:,5),'Uniform',0);
        XX1 = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});
        
        Mat_insurance(jtau1,jtau2)=mean(mean(XX1*Resqtrue_c));
    end
end
% Figure 5c
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);   
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[-0.1 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1))
text(0.8,0.01, 0.85, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 0.75, ['\sigma = ' num2str(std(Mat_insurance(:)),2)])
view([140 15]);
caxis([0 1])
print(figs(it),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_eta_xi_%d_tau%d.eps',model,it),'-depsc');

end

%%% (b) assets derivative wrt eta 
sample=find(D_t);
Mateta_t = Mateta_true(:,1:T-1);
Mateta_t = Mateta_t(sample);
Matc_t = Matc_true(:,1:T-1);
Matc_t=Matc_t(sample);
for it = 1:5
tau_use = quant_plot_grid(it);
Matxi_use = quantile(Matxi_true,tau_use);

Mat_insurance=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        
        
        i_Vect_AGE_t=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
        i_Vect_A_lag=arrayfun(@(x) hermite(x,(Vec1(jtau1)-meanA)/stdA),ua(:,1),'Uniform',0);
        i_Vect_Matdraw_lag=arrayfun(@(x) x./stdY.*hermite(x-1,(Mateta_t-meanY)/stdY),ua(:,2),'Uniform',0);
        i_Vect_Y_lag=arrayfun(@(x) hermite(x,(Y_t-Mateta_t - meanY)/stdY),ua(:,3),'Uniform',0);
        i_Vect_C_lag=arrayfun(@(x) hermite(x,(Matc_t-meanC)/stdC),ua(:,5),'Uniform',0);
        i_Vect_Xi_lag=arrayfun(@(x) hermite(x,(Matxi_use-meanC)/stdC),ua(:,6),'Uniform',0);
        
        i_MatA = cat(2,i_Vect_A_lag{:}).*cat(2,i_Vect_C_lag{:}).*cat(2,i_Vect_Y_lag{:}).*cat(2,i_Vect_Matdraw_lag{:}).*cat(2,i_Vect_AGE_t{:}).*cat(2,i_Vect_Xi_lag{:});
        
        Mat_insurance(jtau1,jtau2)=mean(mean(i_MatA*Resqtrue_a));
    end
end
% Figure 5c
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);   
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Asset response (asset equation)','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[-0.5 2])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(-0.5:0.5:2))
text(0.8,0.01, 1.8, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 1.55, ['\sigma = ' num2str(std(Mat_insurance(:)),2)])
view([140 15]);
caxis([0 1])

print(figs(it),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_a_eta_%d_tau%d.eps',model,it),'-depsc');

end

% (c) asset derivative wrt lag asset 
for it = 1:5
tau_use = quant_plot_grid(it);
Matxi_use = quantile(Matxi_true,tau_use);

Mat_insurance=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        
        
        i_Vect_AGE_t=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
        i_Vect_A_lag=arrayfun(@(x) x./stdA.*hermite(x-1,(Vec1(jtau1)-meanA)/stdA),ua(:,1),'Uniform',0);
        i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(Mateta_t-meanY)/stdY),ua(:,2),'Uniform',0);
        i_Vect_Y_lag=arrayfun(@(x) hermite(x,(Y_t-Mateta_t - meanY)/stdY),ua(:,3),'Uniform',0);
        i_Vect_C_lag=arrayfun(@(x) hermite(x,(Matc_t-meanC)/stdC),ua(:,5),'Uniform',0);
        i_Vect_Xi_lag=arrayfun(@(x) hermite(x,(Matxi_use-meanC)/stdC),ua(:,6),'Uniform',0);
        
        i_MatA = cat(2,i_Vect_A_lag{:}).*cat(2,i_Vect_C_lag{:}).*cat(2,i_Vect_Y_lag{:}).*cat(2,i_Vect_Matdraw_lag{:}).*cat(2,i_Vect_AGE_t{:}).*cat(2,i_Vect_Xi_lag{:});
        
        Mat_insurance(jtau1,jtau2)=mean(mean(i_MatA*Resqtrue_a));
    end
end
% Figure 5c
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);   
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Asset response (asset equation)','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1))
text(0.8,0.01, 0.25, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 0.15, ['\sigma = ' num2str(std(Mat_insurance(:)),2)])
view([140 15]);
caxis([0 1])

print(figs(it),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_a_a_%d_tau%d.eps',model,it),'-depsc');

end


% (ii) PR
Mat_insurance=zeros(Ntau,Ntau);
sample=find(D_t);
Matxi_t= Matxi_true;
for tt=3:T
    Matxi_t=[Matxi_t; Matxi_true];
end
Matxi_t=Matxi_t(sample);

A_lag = Mata_true(:,1:T-1);
A_lag=A_lag(sample);
A_t= Mata_true(:,2:T);
A_t=A_t(sample);
Mateta_lag = Mateta_true(:,1:T-1);
Mateta_lag =Mateta_lag(sample);

[uc1,vc,wc,xc,yc,zc] = ndgrid(0:3);
ua2=[uc1(:),vc(:),wc(:),xc(:),yc(:),zc(:)];
ua2= ua2(ua2(:,1)<=M6,:);
ua2= ua2(ua2(:,2)<=M7,:);
ua2= ua2(ua2(:,3)<=M8,:);
ua2= ua2(ua2(:,4)<=M9,:);
ua2= ua2(ua2(:,5)<=0,:);
ua2= ua2(ua2(:,6)<=M11,:);

Vect_AGE_t=arrayfun(@(x) hermite(x,(AGE_t-meanAGE)/stdAGE),ua2(:,4),'Uniform',0);
Vect_A_lag=arrayfun(@(x) hermite(x,(A_lag-meanA)/stdA),ua2(:,1),'Uniform',0);
Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(Mateta_lag-meanY)/stdY),ua2(:,2),'Uniform',0);
Vect_Xi_lag=arrayfun(@(x) hermite(x,(Matxi_t-meanC)/stdC),ua2(:,6),'Uniform',0);
Vect_Y_lag=arrayfun(@(x) hermite(x,(Y_lag-Mateta_lag-meanY)/stdY),ua2(:,3),'Uniform',0);
MatAlag = cat(2,Vect_A_lag{:}).*cat(2,Vect_Y_lag{:}).*cat(2,Vect_Matdraw_lag{:}).*cat(2,Vect_AGE_t{:}).*cat(2,Vect_Xi_lag{:});

Resqdata_a=zeros(size(MatAlag,2),Ntau);
for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(MatAlag,A_t,tau);
    Resqdata_a(:,jtau)=beta1;
end

for it = 1:5
tau_use = quant_plot_grid(it);
Matxi_use = quantile(Matxi_true,tau_use);

Mat_insurance=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        
        
        Vect_AGE_t=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),ua2(:,4),'Uniform',0);
        Vect_A_lag=arrayfun(@(x) hermite(x,(Vec1(jtau1)-meanA)/stdA),ua2(:,1),'Uniform',0);
        Vect_Matdraw_lag=arrayfun(@(x) x./stdY.*hermite(x-1,(Mateta_lag-meanY)/stdY),ua2(:,2),'Uniform',0);
        Vect_Xi_lag=arrayfun(@(x) hermite(x,(Matxi_use-meanC)/stdC),ua2(:,6),'Uniform',0);
        Vect_Y_lag=arrayfun(@(x) hermite(x,(Y_lag-Mateta_lag-meanY)/stdY),ua2(:,3),'Uniform',0);
        
        i_MatA = cat(2,Vect_A_lag{:}).*cat(2,Vect_Y_lag{:}).*cat(2,Vect_Matdraw_lag{:}).*cat(2,Vect_AGE_t{:}).*cat(2,Vect_Xi_lag{:});
        
        Mat_insurance(jtau1,jtau2)=mean(mean(i_MatA*Resqdata_a));
    end
end
% Figure 5c
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);   
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Asset response (policy)','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[-0.5 1.5])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(-0.5:0.25:1.5))
text(0.8,0.01, 1.25, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 1.1, ['\sigma = ' num2str(std(Mat_insurance(:)),2)])
view([140 15]);
caxis([0 1])

print(figs(it),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_a_eta_pr_%d_tau%d.eps',model,it),'-depsc');

end



for it = 1:5
tau_use = quant_plot_grid(it);
Matxi_use = quantile(Matxi_true,tau_use);

Mat_insurance=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        
        
        Vect_AGE_t=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),ua2(:,4),'Uniform',0);
        Vect_A_lag=arrayfun(@(x) x./stdA.*hermite(x-1,(Vec1(jtau1)-meanA)/stdA),ua2(:,1),'Uniform',0);
        Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(Mateta_lag-meanY)/stdY),ua2(:,2),'Uniform',0);
        Vect_Xi_lag=arrayfun(@(x) hermite(x,(Matxi_use-meanC)/stdC),ua2(:,6),'Uniform',0);
        Vect_Y_lag=arrayfun(@(x) hermite(x,(Y_lag-Mateta_lag-meanY)/stdY),ua2(:,3),'Uniform',0);
        
        i_MatA = cat(2,Vect_A_lag{:}).*cat(2,Vect_Y_lag{:}).*cat(2,Vect_Matdraw_lag{:}).*cat(2,Vect_AGE_t{:}).*cat(2,Vect_Xi_lag{:});
        
        Mat_insurance(jtau1,jtau2)=mean(mean(i_MatA*Resqdata_a));
    end
end
% Figure 5c
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);   
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Asset response (policy)','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1))
text(0.8,0.01, 0.25, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 0.15, ['\sigma = ' num2str(std(Mat_insurance(:)),2)])
view([140 15]);
caxis([0 1])

print(figs(it),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_a_a_pr_%d_tau%d.eps',model,it),'-depsc');

end











sample=find(D(:));
Matxi_tot=Matxi_true;
for tt=2:T
Matxi_tot=[Matxi_tot;Matxi_true];
end
Matxi_tot=Matxi_tot(sample);
Vec1=quantile(Mateta_tot,Vectau);
Vec2=quantile(AGE_tot,Vectau);
Vec3=quantile(Matxi_tot,Vectau);
A_tot = Mata_true(sample);
Mateta_tot =Mateta_true(sample);
Mateps_tot=Mateps_true(sample);


Mat_insurance1=zeros(Ntau,1);
Mat_insurance2=zeros(Ntau,1);
Mat_insurance3=zeros(Ntau,1);

for jtau1=1:Ntau
        Vect_AGE_tot=arrayfun(@(x) hermite(x,(AGE_tot-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_tot=arrayfun(@(x) hermite(x,(A_tot-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_tot=arrayfun(@(x) x*hermite(x-1,(Vec1(jtau1)-meanY)/stdY)./stdY,uc(:,2),'Uniform',0);
        Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_tot-meanY)/stdY),uc(:,3),'Uniform',0);
        Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matxi_tot-meanC)/stdC),uc(:,5),'Uniform',0);
        XX1 = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});
        Mat_insurance1(jtau1)=mean(mean(XX1*Resqtrue_c));
end
for jtau1=1:Ntau
        Vect_AGE_tot=arrayfun(@(x) hermite(x,(Vec2(jtau1)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_tot=arrayfun(@(x) hermite(x,(A_tot-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_tot=arrayfun(@(x) x*hermite(x-1,(Mateta_tot-meanY)/stdY)./stdY,uc(:,2),'Uniform',0);
        Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_tot-meanY)/stdY),uc(:,3),'Uniform',0);
        Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matxi_tot-meanC)/stdC),uc(:,5),'Uniform',0);
        XX1 = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});
        Mat_insurance2(jtau1)=mean(mean(XX1*Resqtrue_c));
end
for jtau1=1:Ntau
        Vect_AGE_tot=arrayfun(@(x) hermite(x,(AGE_tot-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_tot=arrayfun(@(x) hermite(x,(A_tot-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_tot=arrayfun(@(x) x*hermite(x-1,(Mateta_tot-meanY)/stdY)./stdY,uc(:,2),'Uniform',0);
        Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_tot-meanY)/stdY),uc(:,3),'Uniform',0);
        Vect_Xi_tot=arrayfun(@(x) hermite(x,(Vec3(jtau1)-meanC)/stdC),uc(:,5),'Uniform',0);
        XX1 = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});
        Mat_insurance3(jtau1)=mean(mean(XX1*Resqtrue_c));
end


% Figure 5c
figs(1)=figure;
set(figs(1), 'Position', [10 10 500 400]);   
plot(Vectau,Mat_insurance1)
hold on
plot(Vectau,Mat_insurance2)
plot(Vectau,Mat_insurance3)
hold off
xlabel('percentile \tau_{x}','FontSize',8)
ylabel('Average persistence','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
legend('\eta','Age','\xi')
print(figs(1),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_avg_persistence_%d.eps',model),'-depsc');

sample=find(D);
XX=[];
for kk1=0:2
    for kk2=0:2
        XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
            .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
    end    
end

dim=size(XX)
Resqinit_cons_plot=zeros(dim(2),2);
Resqinit_cons_plot(:,1)=rq(XX,Matc_true(sample),0.1);
Resqinit_cons_plot(:,2)=rq(XX,Matc_true(sample),0.9);
figs(1)=figure;
set(figs(1), 'Position', [10 10 500 400]);   
for it = 1:5

tau_use = quant_plot_grid(it);
Matxi_use = quantile(Matxi_true,tau_use);


XX=[];
for kk1=0:2
    for kk2=0:2
        XX=[XX hermite(kk1,(quantile(AGE_tot,Vectau)-meanAGE)/stdAGE)...
            .*hermite(kk2,(Matxi_use-meanC)/stdC)];
    end    
end
iqr = XX*Resqinit_cons_plot(:,2) - XX*Resqinit_cons_plot(:,1)
% Figure 5c

plot(quantile(AGE_tot,Vectau),iqr)
hold on
xlabel('percentile \tau_{age}','FontSize',8)
ylabel('Consumption IQR','FontSize',8)
set(gca,'xlim',[25 60])
set(gca,'ylim',[0.5 0.75])
% set(gca,'xtick',(0:0.2:1))
%set(gca,'ytick',(0:0.2:1))

end
legend({'\tau_{0.1}','\tau_{0.25}','\tau_{0.5}','\tau_{0.75}','\tau_{0.9}'},'Location','SouthWest')
print(figs(1),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_c_iqr_%d.eps',model),'-depsc');


sample=find(D(:));
Matxi_tot=Matxi_true;
for tt=2:T
Matxi_tot=[Matxi_tot;Matxi_true];
end
Matxi_tot=Matxi_tot(sample);
Vect_AGE_tot=arrayfun(@(x) hermite(x,(AGE(sample)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
Vect_A_tot=arrayfun(@(x) hermite(x,(Mata_true(sample)-meanA)/stdA),uc(:,1),'Uniform',0);
Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(Mateta_true(sample)-meanY)/stdY),uc(:,2),'Uniform',0);
Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_true(sample)-meanY)/stdY),uc(:,3),'Uniform',0);
Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matxi_tot-meanC)/stdC),uc(:,5),'Uniform',0);
MatClag = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});

c_bar = MatClag*mean(Resqfinal_cons,2);

sample=find(D);
XX=[];
for kk1=0:2
    for kk2=0:2
        XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
            .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
    end    
end

dim=size(XX)
Resqinit_cons_plot=zeros(dim(2),2);
Resqinit_cons_plot(:,1)=rq(XX,c_bar,0.1);
Resqinit_cons_plot(:,2)=rq(XX,c_bar,0.9);
figs(2)=figure;
set(figs(2), 'Position', [10 10 500 400]);   
for it = 1:5

tau_use = quant_plot_grid(it);
Matxi_use = quantile(Matxi_true,tau_use);


XX=[];
for kk1=0:2
    for kk2=0:2
        XX=[XX hermite(kk1,(quantile(AGE_tot,Vectau)-meanAGE)/stdAGE)...
            .*hermite(kk2,(Matxi_use-meanC)/stdC)];
    end    
end
iqr = XX*Resqinit_cons_plot(:,2) - XX*Resqinit_cons_plot(:,1)
% Figure 5c

plot(quantile(AGE_tot,Vectau),iqr)
hold on
xlabel('percentile \tau_{age}','FontSize',8)
ylabel('Consumption IQR','FontSize',8)
set(gca,'xlim',[25 60])
set(gca,'ylim',[0.15 0.4])
% set(gca,'xtick',(0:0.2:1))
%set(gca,'ytick',(0:0.2:1))

end
legend({'\tau_{0.1}','\tau_{0.25}','\tau_{0.5}','\tau_{0.75}','\tau_{0.9}'},'Location','SouthWest')
print(figs(2),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_c_bar_iqr_%d.eps',model),'-depsc');


% Lifecycle plots
% Figure 5c
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

eta_tot=Mateta_true(sample);

XX=[];
for kk1=0:3
    for kk2=0:3
        for kk2t=0:K2t
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
        end
    end
end
beta = pinv(XX)*eta_tot;


age_list = 25:1:60;



profile=zeros(size(age_list,2),Ntau_plot);

for it = 1:Ntau_plot
tau_use = quant_plot_grid(it);


Matxi_use = quantile(Matxi_true,tau_use);
XX=[];
for kk1=0:3
    for kk2=0:3
            XX=[XX hermite(kk1,(age_list'-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_use-meanC)/stdC)];
    end
end


profile(:,it)=XX*beta;
end

XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)];
end

beta1=rq(XX,C_tot,0.1);
beta9=rq(XX,C_tot,0.9);

XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(age_list'-meanAGE)/stdAGE)];
end
fitted1 = XX*beta1;
fitted9 = XX*beta9;

plot(25:1:60,profile)
legend('\tau_{10}','\tau_{25}','\tau_{50}','\tau_{75}','\tau_{90}')
xlabel('Age','FontSize',8)
ylabel('Fitted value of \eta','FontSize',8)
hold on
plot(25:1:60,[fitted1 fitted9],'--','Color','black');
ylim([-0.6 0.6])
print(figs(40),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/n_xi_profile_%d.eps',model),'-depsc');


% Lifecycle plots
% Figure 5c
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

C_tot=Matc_true(sample);

XX=[];
for kk1=0:3
    for kk2=0:3
        for kk2t=0:K2t
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
        end
    end
end
beta = pinv(XX)*C_tot;


age_list = 25:1:60;



profile=zeros(size(age_list,2),Ntau_plot);

for it = 1:Ntau_plot


tau_use = quant_plot_grid(it);


Matxi_use = quantile(Matxi_true,tau_use);XX=[];
for kk1=0:3
    for kk2=0:3
            XX=[XX hermite(kk1,(age_list'-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_use-meanC)/stdC)];
    end
end


profile(:,it)=XX*beta;
end

XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)];
end

beta1=rq(XX,C_tot,0.1);
beta9=rq(XX,C_tot,0.9);

XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(age_list'-meanAGE)/stdAGE)];
end
fitted1 = XX*beta1;
fitted9 = XX*beta9;
plot(25:1:60,profile)
legend('\tau_{10}','\tau_{25}','\tau_{50}','\tau_{75}','\tau_{90}')
xlabel('Age','FontSize',8)
ylabel('Fitted value of consumption','FontSize',8)
hold on
plot(25:1:60,[fitted1 fitted9],'--','Color','black');
ylim([-0.6 0.6])
print(figs(40),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/c_xi_profile_%d.eps',model),'-depsc');

% Lifecycle plots
% Figure 5c
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

A_tot=Mata_true(sample);

XX=[];
for kk1=0:3
    for kk2=0:3
        for kk2t=0:K2t
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
        end
    end
end
beta = pinv(XX)*A_tot;


age_list = 25:1:60;



profile=zeros(size(age_list,2),Ntau_plot);

for it = 1:Ntau_plot


tau_use = quant_plot_grid(it);


Matxi_use = quantile(Matxi_true,tau_use);

XX=[];
for kk1=0:3
    for kk2=0:3
            XX=[XX hermite(kk1,(age_list'-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_use-meanC)/stdC)];
    end
end


profile(:,it)=XX*beta;
end

XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)];
end

beta1=rq(XX,A_tot,0.1);
beta9=rq(XX,A_tot,0.9);

XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(age_list'-meanAGE)/stdAGE)];
end
fitted1 = XX*beta1;
fitted9 = XX*beta9;
plot(25:1:60,profile)
legend('\tau_{10}','\tau_{25}','\tau_{50}','\tau_{75}','\tau_{90}')
xlabel('Age','FontSize',8)
ylabel('Fitted value of assets','FontSize',8)
hold on
plot(25:1:60,[fitted1 fitted9],'--','Color','black');
print(figs(40),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/a_xi_profile_%d.eps',model),'-depsc');


