% ABBL (2021)
% Code to generate the figures for the baseline model
% NB: assets modelled as function of lagged income only in this baseline

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
model = 11;

%%
% 1. Load results and extra data for plotting
%%%%%%%%%%%%

load('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Quantile model/Results/REVISION_FINAL.mat')

% Check marginal likelihood
figs(1)=figure;
plot(mat_lik)
xlim([0, maxiter])
ylabel('Likelihood','FontSize',9)
xlabel('Iteration','FontSize',9)

% Parent assets + asset residuals
PA=data(:,25);
MatPA=zeros(N,T);
for tt=1:T
    MatPA(:,tt)=PA(tt:T:N*T);
end
PA = MatPA;

PAu=data(:,28);
MatPAu=zeros(N,T);
for tt=1:T
    MatPAu(:,tt)=PAu(tt:T:N*T);
end
PAu = MatPAu;

% Parent income + residuals
PI=data(:,26);
MatPI=zeros(N,T);
for tt=1:T
    MatPI(:,tt)=PI(tt:T:N*T);
end
PI =MatPI;

PIu=data(:,29);
MatPIu=zeros(N,T);
for tt=1:T
    MatPIu(:,tt)=PIu(tt:T:N*T);
end
PIu = MatPIu;

% Parent consumption + residuals
PC=data(:,27);
MatPC=zeros(N,T);
for tt=1:T
    MatPC(:,tt)=PC(tt:T:N*T);
end
PC = MatPC;

PCu=data(:,30);
MatPCu=zeros(N,T);
for tt=1:T
    MatPCu(:,tt)=PCu(tt:T:N*T);
end
PCu = MatPCu;

% Fitted value income
Yhat=data(:,13);
MatY=zeros(N,T);
for tt=1:T
    MatY(:,tt)=Yhat(tt:T:N*T);
end
Yhat=MatY;

% Fitted value consumption
Chat=data(:,15);
MatC=zeros(N,T);
for tt=1:T
    MatC(:,tt)=Chat(tt:T:N*T);
end
Chat=MatC;

% Fitted value assets
Ahat=data(:,18);
MatA=zeros(N,T);
for tt=1:T
    MatA(:,tt)=Ahat(tt:T:N*T);
end
Ahat=MatA;

Nis_scale = 250; % Numerical parameter for the posterior simulations
rng('shuffle'); % Fix seed

% Load the consumption data
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
    
    a_i_old=[]; % Initialize heterogeneity
    Matdraw_old=[]; % Initialize eta draws
    acc=zeros(1,draws); % Stores acceptances
    
    ESS_i=zeros(draws,T); % Stores ESS for person i
    lik_i=zeros(draws,1); % Stores marginal likelihoods for person i
    
    T1=ENT(iii); % Starting time period
    TT=EXT(iii); % Ending time period
    
    Nis=Nis_scale*(TT-T1+1); % No. of importance samples in E-step
    tol=Nis/tol_dnom; % Resampling tolerance
    
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
        
        % Initialize the first MH step
        if j==1
            if iter==1
                a_i_new =i_MatXlag(1,:)*OLSnew_xi + mvnrnd(0,MH_scale*vxi,1);
            elseif iter>1
                a_i_new = xi_draws(iii)
            end
            
            % Update the subsequent MH steps
        elseif j>1
            a_i_new = a_i_old + mvnrnd(0,MH_scale*vxi,1);
        end
        
        % K3t allows for time effects in epsilon (not used in baseline)
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
                
                % Closed form proposal distributions
                alpha=veta1/(veta1 + veps);
                Mvar=2/((1/veta1)+(1/veps));
                Mmu=(1-alpha)*(i_MatAGE1*OLSnew_e0) + alpha*(i_Y(1,ttt));
                
                % Generate the importance samples and their weights
                particle_draw(:,ttt) = Mmu + mvnrnd(0,Mvar,Nis);
                
                % Generate regressors - CONSUMPTION function
                i_Vect_AGE_tot=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
                i_Vect_A_tot=arrayfun(@(x) hermite(x,(i_A(:,ttt)-meanA)/stdA),uc(:,1),'Uniform',0);
                i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
                i_Vect_Eps_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-particle_draw(:,ttt)-meanY)/stdY),uc(:,3),'Uniform',0);
                i_Vect_Xi_tot=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),uc(:,5),'Uniform',0);
                i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:}).*cat(2,i_Vect_Xi_tot{:});
                
                % Generate regressors - ASSETS function
                % (NB: in baseline asset rule is irrelevant and will cancel from likelihood)
                i_Vect_AGE1=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua1(:,2),'Uniform',0);
                i_Vect_Matdraw1=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),ua1(:,1),'Uniform',0);
                i_Vect_Xi1=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),ua1(:,3),'Uniform',0);
                i_Vect_YB=arrayfun(@(x) hermite(x,(YB(iii)-meanYB)/stdYB),ua1(:,4),'Uniform',0);
                i_Vect_ED=arrayfun(@(x) hermite(x,(EDUC(iii)-meanEDUC)/stdEDUC),ua1(:,5),'Uniform',0);
                i_MatA =cat(2,i_Vect_Matdraw1{:}).*cat(2,i_Vect_AGE1{:}).*cat(2,i_Vect_Xi1{:}).*cat(2,i_Vect_YB{:}).*cat(2,i_Vect_ED{:});
                
                % Importance sampling weights
                weights=IS_weight_c_a(i_Y,i_C,particle_draw,i_MatAGE_tot,i_MatAGE1,ttt,...
                    Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,Resqinit_cons,...
                    b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,b1_c,bL_c,...
                    Nis,T1,[],i_MatC,Resqinit_a,b1_a,bL_a,Resqinit_a1,...
                    b1_a1,bL_a1,i_MatA,i_A,asset_rule);
                
                weights=weights./mvnpdf(particle_draw(:,ttt)-Mmu,0,Mvar);
                %                 weights(isnan(weights))=0;
                
                % Marginal likelihood update
                lik=mean(weights);
                
            elseif ttt>T1
                
                % Closed form proposal distributions
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
                
                % Generate regressors - CONSUMPTION function
                i_Vect_AGE_tot=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
                i_Vect_A_tot=arrayfun(@(x) hermite(x,(i_A(:,ttt)-meanA)/stdA),uc(:,1),'Uniform',0);
                i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(particle_draw(:,ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
                i_Vect_Eps_tot=arrayfun(@(x) hermite(x,(i_Y(:,ttt)-particle_draw(:,ttt)-meanY)/stdY),uc(:,3),'Uniform',0);
                i_Vect_Xi_tot=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),uc(:,5),'Uniform',0);
                i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:}).*cat(2,i_Vect_Xi_tot{:});
                
                % Generate regressors - ASSETS function
                i_Vect_AGE_t=arrayfun(@(x) hermite(x,(i_AGE(:,ttt)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
                i_Vect_A_lag=arrayfun(@(x) hermite(x,(i_A(:,ttt-1)-meanA)/stdA),ua(:,1),'Uniform',0);
                i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(i_Y(:,ttt-1)-meanY)/stdY),ua(:,2),'Uniform',0);
                i_Vect_Y_lag=arrayfun(@(x) hermite(x,(i_Y(:,ttt-1)-particle_draw(:,ttt-1) - meanY)/stdY),ua(:,3),'Uniform',0);
                i_Vect_C_lag=arrayfun(@(x) hermite(x,(i_C(:,ttt-1)-meanC)/stdC),ua(:,5),'Uniform',0);
                i_Vect_Xi_lag=arrayfun(@(x) hermite(x,(a_i_new-meanC)/stdC),ua(:,6),'Uniform',0);
                i_MatA = cat(2,i_Vect_A_lag{:}).*cat(2,i_Vect_C_lag{:}).*cat(2,i_Vect_Y_lag{:}).*cat(2,i_Vect_Matdraw_lag{:}).*cat(2,i_Vect_AGE_t{:}).*cat(2,i_Vect_Xi_lag{:});
                
                % Importance sampling weights
                weights=weights...
                    .*IS_weight_c_a(i_Y,i_C,particle_draw,i_MatAGE_tot,i_MatAGE1,ttt,...
                    Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,Resqinit_cons,...
                    b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,b1_c,bL_c,...
                    Nis,T1,i_Matdraw_lag,i_MatC,Resqinit_a,b1_a,bL_a,Resqinit_a1,b1_a1,bL_a1,i_MatA,i_A,asset_rule);
                
                weights=weights./mvnpdf(particle_draw(:,ttt)-Mmu,0,Mvar);
                % weights(isnan(weights))=0;
                
                % Marginal likelihood update
                lik=lik*sum(weights);
                
            end
            
            % Self-normalization
            weights=weights./sum(weights);
            
            % Effective sample size
            ESS= 1/sum(weights.^2);
            ESS_i(j,ttt)=ESS;
            
            % Adaptive resampling (tolerance set by tol)
            if ttt<TT && ESS<tol && ~isnan(ESS)
                particle_draw = datasample(particle_draw,Nis,'weights',weights);
                weights=ones(Nis,1)*(1/Nis);
            elseif ttt==TT && ~isnan(ESS)
                Matdraw_new=datasample(particle_draw,1,'weights',weights);
                
                % (VERY RARE)
                % hacky way to reset things if we get an
                % underflow issue - use the draw from last MH iter
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
        
        % PMCMC acceptance object
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
    
    % Store generated samples and diagnostics
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

%% 
% 3. Motivating figures
%%%%%%%%%%%%

%%% Motivating figures in main paper

Y_bar = [];
C_bar = [];
A_bar =[];
AGE_bar = [];
for iii=1:N
    length = EXT(iii)-ENT(iii);
    Y_bar = [Y_bar;mean(Y(iii,ENT(iii):EXT(iii)))];
    C_bar = [C_bar;mean(C(iii,ENT(iii):EXT(iii)))];
    A_bar = [A_bar; mean(A(iii,ENT(iii):EXT(iii)))];
    AGE_bar = [AGE_bar;mean(AGE(iii,ENT(iii):EXT(iii)))];
end

XXX=[ones(size(Y_bar)) Y_bar AGE_bar A_bar C_bar];
XXX_all =[ones(size(Y_tot)) Y_tot AGE_tot A_tot C_tot];

% Non-linear - averages
K=2;
X=[];
for kk1=0:K
    for kk2=0:K
        for kk3=0:(K-1)
            X=[X hermite(kk1,(Y_bar-meanY)/stdY)...
                .*hermite(kk2,(A_bar-meanA)/stdA)...
                .*hermite(kk3,(AGE_bar-meanAGE)/stdAGE)];
        end
    end
end

Ntau2=95;
Vectau2=linspace(0,1,Ntau2);
Vectau2=Vectau2(3:end-3);

Resqdata_cons=zeros(size(X,2),size(Vectau2,2));
for jtau=1:size(Vectau2,2)
    beta1=rq(X,C_bar,Vectau2(jtau));
    Resqdata_cons(:,jtau)=beta1;
end

Vec1=quantile(A_bar,Vectau);
Vec2=quantile(AGE_bar,Vectau);

Mat_insurance1=zeros(Ntau,Ntau);
Mat_insurance2=zeros(Ntau,Ntau);
Mat_insurance3=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        X=[];
        for kk1=0:K
            for kk2=0:K
                for kk3=0:(K-1)
                    if kk1 <1
                        X=[X zeros(size(Y_bar,1),1)];
                    else
                        X=[X kk1/stdY*hermite(kk1-1,(Y_bar-meanY)/stdY)...
                            .*hermite(kk2,(Vec1(jtau1)-meanA)/stdA)...
                            .*hermite(kk3,(Vec2(jtau2)-meanAGE)/stdAGE)];
                    end
                end
            end
        end     
        Mat_insurance1(jtau1,jtau2)=mean(mean(X*Resqdata_cons(:,1:30)));
        Mat_insurance2(jtau1,jtau2)=mean(mean(X*Resqdata_cons(:,1:90)));
        Mat_insurance3(jtau1,jtau2)=mean(mean(X*Resqdata_cons(:,61:90)));       
    end
end

% Figure 1
figs(4)=figure;
set(figs(4), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_insurance2)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0.2 0.8])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0.2:0.2:8))
ytickangle(-0)
caxis([0 1])
text(0.8,0.01, 0.75, ['\mu = ' num2str(mean(Mat_insurance2(:)),2)])
text(0.8,0.01, 0.675, ['\sigma = ' num2str(std(Mat_insurance2(:)),1)])
view([140 15]);
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/motivate_%d_q2.eps',model),'-depsc');


% Figure 2a
figs(4)=figure;
set(figs(4), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_insurance1)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0.2 0.8])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0.2:0.2:8))
ytickangle(0)
text(0.8,0.01, 0.75, ['\mu = ' num2str(mean(Mat_insurance1(:)),2)])
text(0.8,0.01, 0.675, ['\sigma = ' num2str(std(Mat_insurance1(:)),1)])
view([140 15]);
caxis([0 1])
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/motivate_%d_q1.eps',model),'-depsc');

% Figure 2b
figs(4)=figure;
set(figs(4), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_insurance3)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0.2 0.8])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0.2:0.2:8))
ytickangle(-0)
text(0.8,0.01, 0.75, ['\mu = ' num2str(mean(Mat_insurance3(:)),2)])
text(0.8,0.01, 0.675, ['\sigma = ' num2str(std(Mat_insurance3(:)),1)])
view([140 15]);
caxis([0 1])
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/motivate_%d_q3.eps',model),'-depsc');

%%% Bootstrap figures for appendix

XXX=[ones(size(Y_bar)) Y_bar AGE_bar A_bar C_bar];

Msim=250;

Mat_insurance1=zeros(Ntau,Ntau,Msim);
Mat_insurance2=zeros(Ntau,Ntau,Msim);
Mat_insurance3=zeros(Ntau,Ntau,Msim);

Vec1=quantile(A_bar,Vectau);
Vec2=quantile(AGE_bar,Vectau);

for nboot =1:Msim
nboot    
XX_boot = datasample(XXX,N);
Y_bar = XX_boot(:,2);
AGE_bar = XX_boot(:,3);
A_bar = XX_boot(:,4);
C_bar = XX_boot(:,5);

K=2;
X=[];
for kk1=0:K
    for kk2=0:K
        for kk3=0:(K-1)
            X=[X hermite(kk1,(Y_bar-meanY)/stdY)...
                .*hermite(kk2,(A_bar-meanA)/stdA)...
                .*hermite(kk3,(AGE_bar-meanAGE)/stdAGE)];
        end
    end
end

Ntau2=47;
Vectau2=linspace(0,1,Ntau2);
Vectau2=Vectau2(3:end-3);

Resqdata_cons=zeros(size(X,2),size(Vectau2,2));
for jtau=1:size(Vectau2,2)
    beta1=rq(X,C_bar,Vectau2(jtau));
    Resqdata_cons(:,jtau)=beta1;
end

for jtau1=1:Ntau
    for jtau2=1:Ntau       
        X=[];
        for kk1=0:K
            for kk2=0:K
                for kk3=0:(K-1)
                    if kk1 <1
                        X=[X zeros(size(XXX,1),1)];
                    else
                        X=[X kk1/stdY*hermite(kk1-1,(Y_bar-meanY)/stdY)...
                            .*hermite(kk2,(Vec1(jtau1)-meanA)/stdA)...
                            .*hermite(kk3,(Vec2(jtau2)-meanAGE)/stdAGE)];
                    end
                end
            end
        end       
        Mat_insurance1(jtau1,jtau2,nboot)=mean(mean(X*Resqdata_cons(:,1:14)));
        Mat_insurance2(jtau1,jtau2,nboot)=mean(mean(X*Resqdata_cons(:,1:42)));
        Mat_insurance3(jtau1,jtau2,nboot)=mean(mean(X*Resqdata_cons(:,29:42)));
    end
end
end

% Figure A3a
figs(4)=figure;
Mat1 = quantile(Mat_insurance2,0.1,3);
Mat2 = quantile(Mat_insurance2,0.9,3);
surf(Vectau,Vectau,Mat1,'FaceAlpha',0.4)
hold on
surf(Vectau,Vectau,Mat2,'FaceAlpha',0.2)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.25:1))
view([140 15]);
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/motivate_boot_q2.eps'),'-depsc');

% Figure A3b
figs(4)=figure;
Mat1 = quantile(Mat_insurance1,0.1,3);
Mat2 = quantile(Mat_insurance1,0.9,3);
surf(Vectau,Vectau,Mat1,'FaceAlpha',0.4)
hold on
surf(Vectau,Vectau,Mat2,'FaceAlpha',0.2)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.25:1))
view([140 15]);
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/motivate_boot_q1.eps'),'-depsc');

% Figure A3c
figs(4)=figure;
Mat1 = quantile(Mat_insurance3,0.1,3);
Mat2 = quantile(Mat_insurance3,0.9,3);
surf(Vectau,Vectau,Mat1,'FaceAlpha',0.4)
hold on
surf(Vectau,Vectau,Mat2,'FaceAlpha',0.2)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.25:1))
view([140 15]);
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/motivate_boot_q3.eps'),'-depsc');

%%
% 4. Main figures 
%%%%%%%%%%%%

%%% FIGURE 4a - Main Insurance Y
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
print(figs(4),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/figs/abb_cy_single_%d.eps',model),'-depsc');

%%% FIGURE 4d - Average Insurance

sample=find(D(:));
Vec1=quantile(A_tot(:),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

Mat_insurance=zeros(Ntau,Ntau);

Matxi_tot= Matdraw(:,T+1);
for tt=2:T
    Matxi_tot=[Matxi_tot; Matdraw(:,T+1)];
end
sample=find(D);
Mateta_tot = Matdraw(:,1:T);
Mateta_tot =Mateta_tot(sample);
Mateps_tot = Y_tot - Mateta_tot;
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
print(figs(4),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/abb_c_single_%d.eps',model),'-depsc');


%%
% 5. Main heterogeneity results
%%%%%%%%%%%%

%%% FIGURE 5 - consumption response by types

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
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);
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
print(figs(it),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/abb_eta_xi_%d_tau%d.eps',model,it),'-depsc');
end

%%% FIGURE A19 - robustness to knots

close all
quant_plot_grid = [0.1 0.25 0.5 0.75 0.9];
% pars = load('/home/jdlight/ABBL_PMCMC/JOE_codes/Knot_robustness/Cons/Results/inc_estspec_robust_5.mat')

Vec1=quantile(A_tot(:),Vectau);
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
        Mat_insurance(jtau1,jtau2)=mean(mean(XX1*pars.Resqfinal_cons));
    end
end
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);
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
print(figs(it),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/abb_eta_xi_19_knots_tau%d.eps',it),'-depsc');
end

%%
% 6. Lifecycle results
%%%%%%%%%%%%

%%% CONSUMPTION PROFILES

sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

C_tot=Chat(sample) + Matc_true(sample);

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

%%% Figure 6a
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
plot(25:1:60,profile)
xlabel('Age','FontSize',8)
ylabel('Fitted value of consumption','FontSize',8)
hold on
plot(25:1:60,[fitted1 fitted9],'--','Color','black');
legend('\tau_{10}','\tau_{25}','\tau_{50}','\tau_{75}','\tau_{90}','Location','northeastoutside')
print(figs(40),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/chat_xi_profile_%d.eps',model),'-depsc');

% ASSET PROFILES

sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

A_tot=Ahat(sample) + Mata_true(sample);

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

%%% Figure 6b
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
plot(25:1:60,profile)
xlabel('Age','FontSize',8)
ylabel('Fitted value of assets','FontSize',8)
hold on
plot(25:1:60,[fitted1 fitted9],'--','Color','black');
legend('\tau_{10}','\tau_{25}','\tau_{50}','\tau_{75}','\tau_{90}','Location','northeastoutside')
print(figs(40),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/ahat_xi_profile_%d.eps',model),'-depsc');

% ETA PROFILES

sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

ETA_tot=Yhat(sample) + Mateta_true(sample);

XX=[];
for kk1=0:3
    for kk2=0:3
        for kk2t=0:K2t
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
        end
    end
end
beta = pinv(XX)*ETA_tot;

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

beta1=rq(XX,ETA_tot,0.1);
beta9=rq(XX,ETA_tot,0.9);

XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(age_list'-meanAGE)/stdAGE)];
end
fitted1 = XX*beta1;
fitted9 = XX*beta9;

%%% Figure 6c
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
plot(25:1:60,profile)
xlabel('Age','FontSize',8)
ylabel('Fitted value of \eta','FontSize',8)
hold on
plot(25:1:60,[fitted1 fitted9],'--','Color','black');
legend('\tau_{10}','\tau_{25}','\tau_{50}','\tau_{75}','\tau_{90}','Location','northeastoutside')
print(figs(40),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/nhat_xi_profile_%d.eps',model),'-depsc');

%%% FIGURE A3

Ntau_plot=3;
quant_plot_grid = [0.1 0.5 0.9];

sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

c_tot=Chat(sample) + Matc_true(sample);

XX=[];
for kk1=0:3
    for kk2=0:3
        for kk2t=0:K2t
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
        end
    end
end
beta = pinv(XX)*c_tot;
beta1=rq(XX,c_tot,0.1);
beta9=rq(XX,c_tot,0.9);

age_list = 25:1:60;

% FIGURE A12
var_of_c = var(c_tot)
XX=[];
for kk1=0:3
            XX=[XX hermite(kk1,(Matxi_tot-meanC)/stdC)];
end
b_ols = pinv(XX)*c_tot;
c_hat = XX*b_ols;
var_of_c_hat = var(c_hat)

profile=zeros(size(age_list,2),Ntau_plot);
profile_p1=zeros(size(age_list,2),Ntau_plot);
profile_p3=zeros(size(age_list,2),Ntau_plot);

for it = 1:3
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
profile_p1(:,it) = XX*beta1;
profile_p3(:,it) = XX*beta9;
end

% figure A12
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
subplot(1,3,1)
hold on
plot(25:1:60,profile(:,1),'Color','red')
plot(25:1:60,profile_p1(:,1),'--','Color','red')
plot(25:1:60,profile_p3(:,1),'--','Color','red')
ylim([9 12])
xlabel('Age','FontSize',8)
ylabel('Fitted value of consumption','FontSize',8)
subplot(1,3,2)
hold on
plot(25:1:60,profile(:,2),'Color','blue')
plot(25:1:60,profile_p1(:,2),'--','Color','blue')
plot(25:1:60,profile_p3(:,2),'--','Color','blue')
ylim([9 12])
xlabel('Age','FontSize',8)
ylabel('Fitted value of consumption','FontSize',8)
subplot(1,3,3)
hold on
plot(25:1:60,profile(:,3),'Color','green')
plot(25:1:60,profile_p1(:,3),'--','Color','green')
plot(25:1:60,profile_p3(:,3),'--','Color','green')
ylim([9 12])
xlabel('Age','FontSize',8)
ylabel('Fitted value of consumption','FontSize',8)
print(figs(40),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/chat_xi_profilec_%d.eps',model),'-depsc');

%%% FIGURE A13
Ntau_plot=3;
quant_plot_grid = [0.1 0.5 0.9];

sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

a_tot=Ahat(sample) + Mata_true(sample);

XX=[];
for kk1=0:3
    for kk2=0:3
        for kk2t=0:K2t
            XX=[XX hermite(kk1,(AGE_tot-meanAGE)/stdAGE)...
                .*hermite(kk2,(Matxi_tot-meanC)/stdC)];
        end
    end
end
beta = pinv(XX)*a_tot;
beta1=rq(XX,a_tot,0.1);
beta9=rq(XX,a_tot,0.9);

age_list = 25:1:60;
profile=zeros(size(age_list,2),Ntau_plot);
profile_p1=zeros(size(age_list,2),Ntau_plot);
profile_p3=zeros(size(age_list,2),Ntau_plot);

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
profile_p1(:,it) = XX*beta1;
profile_p3(:,it) = XX*beta9;
end

% figure A13a
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
hold on
plot(25:1:60,profile(:,1),'Color','red')
plot(25:1:60,profile(:,2),'Color','blue')
plot(25:1:60,profile(:,3),'Color','green')
xlabel('Age','FontSize',8)
ylabel('Fitted value of assets','FontSize',8)
plot(25:1:60,profile_p1(:,1),'--','Color','red')
plot(25:1:60,profile_p1(:,2),'--','Color','blue')
plot(25:1:60,profile_p1(:,3),'--','Color','green')
plot(25:1:60,profile_p3(:,1),'--','Color','red')
plot(25:1:60,profile_p3(:,2),'--','Color','blue')
plot(25:1:60,profile_p3(:,3),'--','Color','green')
legend('\tau_{10}','\tau_{50}','\tau_{90}','Location','northeastoutside')
print(figs(40),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/ahat_xi_profileb_%d.eps',model),'-depsc');

%%% FIGURE A13b

Ntau_plot=3;
quant_plot_grid = [0.1 0.5 0.9];

sample=find(D);
Matxi_tot=Matxi_true;
for tt=2:T
   Matxi_tot=[Matxi_tot;Matxi_true]; 
end
Matxi_tot=Matxi_tot(sample);

eta_tot=Yhat(sample) + Mateta_true(sample);

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
beta1=rq(XX,eta_tot,0.1);
beta9=rq(XX,eta_tot,0.9);

age_list = 25:1:60;
profile=zeros(size(age_list,2),Ntau_plot);
profile_p1=zeros(size(age_list,2),Ntau_plot);
profile_p3=zeros(size(age_list,2),Ntau_plot);

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
profile_p1(:,it) = XX*beta1;
profile_p3(:,it) = XX*beta9;
end

% figure A13b
figs(40)=figure;
set(figs(40), 'Position', [10 10 800 300]);
hold on
plot(25:1:60,profile(:,1),'Color','red')
plot(25:1:60,profile(:,2),'Color','blue')
plot(25:1:60,profile(:,3),'Color','green')
xlabel('Age','FontSize',8)
ylabel('Fitted value of \eta','FontSize',8)
plot(25:1:60,profile_p1(:,1),'--','Color','red')
plot(25:1:60,profile_p1(:,2),'--','Color','blue')
plot(25:1:60,profile_p1(:,3),'--','Color','green')
plot(25:1:60,profile_p3(:,1),'--','Color','red')
plot(25:1:60,profile_p3(:,2),'--','Color','blue')
plot(25:1:60,profile_p3(:,3),'--','Color','green')
legend('\tau_{10}','\tau_{50}','\tau_{90}','Location','northeastoutside')
print(figs(40),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/nhat_xi_profileb_%d.eps',model),'-depsc');

%%
% 7. Other Results and Robustness
%%%%%%%%%%%%

%%% FIGURE A11 - stdandard deviation of consumption residuals by types
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
Mat_R2=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau  
        Vect_AGE_tot=arrayfun(@(x) hermite(x,(Vec2(jtau2)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
        Vect_A_tot=arrayfun(@(x) hermite(x,(Vec1(jtau1)-meanA)/stdA),uc(:,1),'Uniform',0);
        Vect_Matdraw_tot=arrayfun(@(x) x*hermite(x-1,(Mateta_tot-meanY)/stdY)./stdY,uc(:,2),'Uniform',0);
        Vect_Eps_tot=arrayfun(@(x) hermite(x,(Mateps_tot-meanY)/stdY),uc(:,3),'Uniform',0);
        Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matxi_use-meanC)/stdC),uc(:,5),'Uniform',0);
        XX1 = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});       
        Insurance_tot = XX1*Resqfinal_cons;
        Insurance_tot_mean = mean(Insurance_tot,2);

        res = Insurance_tot - Insurance_tot_mean;
        
        Mat_insurance(jtau1,jtau2)=std(res(:));
        Mat_R2(jtau1,jtau2) = 1 - var(res(:))/var(Insurance_tot(:));
    end
end


% figure A11
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_R2)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('R2 of consumption response explained by covariates','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1))
text(0.8,0.01, 0.25, ['\mu = ' num2str(mean(Mat_R2(:)),2)])
text(0.8,0.01, 0.15, ['\sigma = ' num2str(std(Mat_R2(:)),1)])
caxis([0 1])
view([140 15]);
print(figs(it),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/r2_eta_xi_%d_tau%d.eps',model,it),'-depsc');
end



%%% FIGURE A25 - QQ plots by type

% Education
figs(13) = figure;
set(figs(13), 'Position', [10 10 600 300]);
sample=find(EDUC(:)==0);
Vect1=Matxi_true(sample);
sample=find(EDUC(:)==1);
Vect2=Matxi_true(sample);
qqplot(Vect1,Vect2)
xlim([-1 1])
ylim([-1 1])
xlabel('Non-graduates')
ylabel('Graduates')
print(figs(13),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/qq_educ_xi_%d.eps',model),'-depsc');

% Cohort 
figs(13) = figure;
set(figs(13), 'Position', [10 10 600 300]);
sample=find(YB(:)<=1969);
Vect1=Matxi_true(sample);
sample=find(YB(:)>=1969);
Vect2=Matxi_true(sample);
qqplot(Vect1,Vect2)
xlim([-1 1])
ylim([-1 1])
xlabel('Older cohorts')
ylabel('Younger cohorts')
print(figs(13),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/qq_coh_xi_%d.eps',model),'-depsc');




%%% FIGURE A26

sample=find(D(:));
Vec1=quantile(Mata_true(sample),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

EDUC_tot=EDUC;
for tt=2:T
    EDUC_tot=[EDUC_tot;EDUC];
end
EDUC_tot=EDUC_tot(sample);

YB_tot=YB;
for tt=2:T
    YB_tot=[YB_tot;YB];
end
YB_tot = YB_tot(sample);

Matxi_tot= Matxi_true;
for tt=2:T
    Matxi_tot=[Matxi_tot; Matxi_true];
end
Matxi_tot = Matxi_tot(sample);

AGE_tot=AGE(sample);
A_tot=A(sample);
Y_tot= Y(sample);
Matdraw_tot = Mateta_true(sample);
C_tot=C(sample);
M6=1;
Resqinit_cons_educ=zeros((M1+1)*(M2+1)*(M3+1)*(M4+1)*(M5+1)*(M6+1),Ntau);
XX=[];

for kk1=0:M1
    for kk2=0:M2
        for kk3=0:M3
            for kk4=0:M4
                for kk5=0:M5
                    for kk6=0:M6
                XX=[XX hermite(kk1,(A_tot-meanA)/stdA)...
                    .*hermite(kk2,(Matdraw_tot-meanY)/stdY)...
                    .*hermite(kk3,(Y_tot-Matdraw_tot-meanY)/stdY)...
                    .*hermite(kk4,(AGE_tot-meanAGE)/stdAGE)...
                    .*hermite(kk5,(Matxi_tot-meanC)/stdC)...
                    .*hermite(kk6,EDUC_tot)];
                    end
                end
            end
        end
    end
end

for jtau=1:Ntau
    tau=Vectau(jtau);
    beta1=rq(XX,C_tot,tau);
    Resqinit_cons_educ(:,jtau)=beta1;
end

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
        
        XX=[];
        for kk1=0:M1
            for kk2=0:M2
                for kk3=0:M3
                    for kk4=0:M4
                        for kk5=0:M5
                            for kk6=0:M6
                                if kk2<1
                                XX = [XX zeros(size(Y_tot))];
                                else
                                XX=[XX hermite(kk1,(Vec1(jtau1)-meanA)/stdA)...
                                    .*kk2.*hermite(kk2-1,(Matdraw_tot-meanY)/stdY)./stdY...
                                    .*hermite(kk3,(Y_tot-Matdraw_tot-meanY)/stdY)...
                                    .*hermite(kk4,(Vec2(jtau2)-meanAGE)/stdAGE)...
                                    .*hermite(kk5,(Matxi_use-meanC)/stdC)...
                                    .*hermite(kk6,EDUC_tot)];
                                end
                            end
                        end
                    end
                end
            end
        end
        
        Mat_insurance(jtau1,jtau2)=mean(mean(XX*Resqinit_cons_educ));
    end
end

% figures by type
figs(it)=figure;
set(figs(it), 'Position', [10 10 500 400]);
surf(Vectau,Vectau,Mat_insurance)
xlabel('Age','FontSize',8)
ylabel('Assets','FontSize',8)
zlabel('Consumption response','FontSize',8)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.25:1))
text(0.8,0.01, 0.85, ['\mu = ' num2str(mean(Mat_insurance(:)),2)])
text(0.8,0.01, 0.75, ['\sigma = ' num2str(std(Mat_insurance(:)),1)])
caxis([0 1])
view([140 15]);
print(figs(it),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/abb_eta_xi_educ_%d_tau%d.eps',model,it),'-depsc');

end

% FIGURE A26d

sample=find(D(:));
Vec1=quantile(A_tot(:),Vectau);
Vec2=quantile(AGE_tot(:),Vectau);

Matxi_tot= Matdraw(:,T+1);
for tt=2:T
    Matxi_tot=[Matxi_tot; Matdraw(:,T+1)];
end
sample=find(D);
Mateta_tot = Matdraw(:,1:T);
Mateta_tot =Mateta_tot(sample);
Mateps_tot = Y_tot - Mateta_tot;
%Mateps_tot=Mateps_true(sample);
Matxi_tot=Matxi_tot(sample);
Mat_insurance=zeros(Ntau,Ntau);

for jtau1=1:Ntau
    for jtau2=1:Ntau
        
        XX=[];
        for kk1=0:M1
            for kk2=0:M2
                for kk3=0:M3
                    for kk4=0:M4
                        for kk5=0:M5
                            for kk6=0:M6
                                if kk2<1
                                XX = [XX zeros(size(Y_tot))];
                                else
                                XX=[XX hermite(kk1,(Vec1(jtau1)-meanA)/stdA)...
                                    .*kk2.*hermite(kk2-1,(Matdraw_tot-meanY)/stdY)./stdY...
                                    .*hermite(kk3,(Y_tot-Matdraw_tot-meanY)/stdY)...
                                    .*hermite(kk4,(Vec2(jtau2)-meanAGE)/stdAGE)...
                                    .*hermite(kk5,(Matxi_tot-meanC)/stdC)...
                                    .*hermite(kk6,EDUC_tot)];
                                end
                            end
                        end
                    end
                end
            end
        end
        
        Mat_insurance(jtau1,jtau2)=mean(mean(XX*Resqinit_cons_educ));
    end
end

% figure 26d
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
print(figs(4),sprintf('/home/jdlight/ABBL_PMCMC/JOE_codes/Replication/Cons/Updated_figs/abb_c_single_educ_%d.eps',model),'-depsc');



%% 
% 8. Motivating regressions with parental variables
%%%%%%%%%%%%

% Expand the sample
Nsim=10;
N=N*Nsim;
AGE=kron(ones(Nsim,1),AGE);
D=kron(ones(Nsim,1),D);
Y=kron(ones(Nsim,1),Y);
C=kron(ones(Nsim,1),C);
A=kron(ones(Nsim,1),A);
EDUC=kron(ones(Nsim,1),EDUC);
YB=kron(ones(Nsim,1),YB);
Chat=kron(ones(Nsim,1),Chat);
Ahat=kron(ones(Nsim,1),Ahat);
MatPC=kron(ones(Nsim,1),MatPC);
MatPA=kron(ones(Nsim,1),MatPA);
MatPI=kron(ones(Nsim,1),MatPI);
MatPCu=kron(ones(Nsim,1),MatPCu);
MatPAu=kron(ones(Nsim,1),MatPAu);
MatPIu=kron(ones(Nsim,1),MatPIu);
xi_draws=kron(ones(Nsim,1),xi_draws);
%%%%

tempdata=zeros(N*T,13);

MatC=zeros(N*T,1);
for ii=1:N
    MatC((ii-1)*T+1:ii*T,:)=C(ii,:);
end
tempdata(:,9)=MatC;

MatA=zeros(N*T,1);
for ii=1:N
    MatA((ii-1)*T+1:ii*T,:)=A(ii,:);
end
tempdata(:,12)=MatA;

MatY=zeros(N*T,1);
for ii=1:N
    MatY((ii-1)*T+1:ii*T,:)=Y(ii,:);
end
tempdata(:,5)=MatY;

MatAGE=zeros(N*T,1);
for ii=1:N
    MatAGE((ii-1)*T+1:ii*T,:)=AGE(ii,:);
end
tempdata(:,3)=MatAGE;

MatD=zeros(N*T,1);
for ii=1:N
    MatD((ii-1)*T+1:ii*T,:)=D(ii,:);
end
tempdata(:,4)=MatD;

MatEDUC=zeros(N*T,1);
for ii=1:N
    MatEDUC((ii-1)*T+1:ii*T,:)=EDUC(ii);
end
tempdata(:,6)=MatEDUC;

MatYB=zeros(N*T,1);
for ii=1:N
    MatYB((ii-1)*T+1:ii*T,:)=YB(ii);
end
tempdata(:,7)=MatYB;

data=tempdata;
clear tempdata
%%%%

load_data_plot(data)
load_cons_data_plot(data)

PA=zeros(N,1);
PI=zeros(N,1);
PC=zeros(N,1);
PAu=zeros(N,1);
PIu=zeros(N,1);
PCu=zeros(N,1);
for iii=1:N
    PAu(iii)=MatPAu(iii,ENT(iii));
    PIu(iii)=MatPIu(iii,ENT(iii));
    PCu(iii)=MatPCu(iii,ENT(iii));
    PA(iii)=MatPA(iii,ENT(iii));
    PI(iii)=MatPI(iii,ENT(iii));
    PC(iii)=MatPC(iii,ENT(iii));
end

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
                i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(i_Y(:,ttt-1)-meanY)/stdY),ua(:,2),'Uniform',0);
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

% Full sample, non-residualized
sample=find((PC~=0).*(PA~=0).*(PI~=0));
XX=[PC(sample) PI(sample) PA(sample)];
g_tot = (1:(N/Nsim))';
for tt=2:Nsim
g_tot=[g_tot;(1:(N/Nsim))'];
end
g_tot=g_tot(sample);
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PC(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PI(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PA(sample), g_tot)

% Young sample, non-residualized
sample=find((AGE1<=45).*(PC~=0).*(PA~=0).*(PI~=0));
XX=[PC(sample) PI(sample) PA(sample)];
g_tot = (1:(N/Nsim))';
for tt=2:Nsim
g_tot=[g_tot;(1:(N/Nsim))'];
end
g_tot=g_tot(sample);
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PC(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PI(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PA(sample), g_tot)

% Full sample, residualized
sample=find((PCu~=0).*(PAu~=0).*(PIu~=0));
XX=[PCu(sample) PIu(sample) PAu(sample)];
g_tot = (1:(N/Nsim))';
for tt=2:Nsim
g_tot=[g_tot;(1:(N/Nsim))'];
end
g_tot=g_tot(sample);
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PCu(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PIu(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PAu(sample), g_tot)

% Full sample, residualized, controlling for YB and EDUC
sample=find((PCu~=0).*(PAu~=0).*(PIu~=0));
XX=[PCu(sample) PIu(sample) PAu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2];
g_tot = (1:(N/Nsim))';
for tt=2:Nsim
g_tot=[g_tot;(1:(N/Nsim))'];
end
g_tot=g_tot(sample);
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)

XX=[PCu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2];
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)

XX=[PIu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2];
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)

XX=[PAu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2]; 
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)

% Young sample, residualized
sample=find((AGE1<=45).*(PCu~=0).*(PAu~=0).*(PIu~=0));
XX=[PCu(sample) PIu(sample) PAu(sample)];
g_tot = (1:(N/Nsim))';
for tt=2:Nsim
g_tot=[g_tot;(1:(N/Nsim))'];
end
g_tot=g_tot(sample);
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PCu(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PIu(sample), g_tot)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), PAu(sample), g_tot)

% Young sample, residualized, controlling for YB and EDUC
sample=find((AGE1<=45).*(PCu~=0).*(PAu~=0).*(PIu~=0));
XX=[PCu(sample) PIu(sample) PAu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2];
g_tot = (1:(N/Nsim))';
for tt=2:Nsim
g_tot=[g_tot;(1:(N/Nsim))'];
end
g_tot=g_tot(sample);
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)
XX=[PCu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2];
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)

XX=[PIu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2];
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)

XX=[PAu(sample) EDUC(sample) AGE1(sample) AGE1(sample).^2]; 
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g_tot)


XX=[ones(size(PC(sample))) PA(sample)];
b_ols=pinv(XX)*Matxi_true(sample)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g(sample,:));
sqrt(diag(varBhat))

XX=[ones(size(PC(sample))) PI(sample)];
b_ols=pinv(XX)*Matxi_true(sample)
[res, varBhat, Rsq] = clusterreg(Matxi_true(sample), XX, g(sample,:));
sqrt(diag(varBhat))
scatter(Matxi_true,PI)





