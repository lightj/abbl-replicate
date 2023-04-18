% ABBL (2023)
% Code to generate the figures for the model with asset rule

cd '/home/jdlight/ABBL_PMCMC/JOE_codes/SMC/Consumption/'

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

% load('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Cons_parametric/ResultsA/REVISION_ASSETS_FINAL.mat')

Resqinit= Resqfinal;
Resqinit_e0 = Resqfinal_e0;
Resqinit_eps=Resqfinal_eps;

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


%%% FIGURE A14d 

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

% Figure A14d
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


%%% FIGURE A14
Ntau_plot=5;
quant_plot_grid = [0.1 0.25 0.5 0.75 0.9];

sample=find(D(:));
Mateta_tot =Mateta_true(sample);
Mateps_tot=Mateps_true(sample);

Vec1=quantile(Mata_true(sample),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

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
% print(figs(it),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_eta_xi_%d_tau%d.eps',model,it),'-depsc');

end


%%% FIGURE A16
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

%%% FIGURE A15
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

