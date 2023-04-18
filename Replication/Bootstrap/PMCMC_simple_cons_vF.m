function [] = PMCMC_simple_cons_FIXED_vF(data,maxitercons,savepathincome,savepathconsumption)
% ABBL (2023)
% Code to estimate the baseline non-linear consumption model
% NB: assets modelled as function of lagged income only in this baseline

%clear all
%clc;

cd ('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Cons/Quantile model')


global T N K1 K2 K2t K3 K4 K5 K6 Ntau Vectau tau AGE1 Y1 A1...
    D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
     Matdraw Matdraw1 MatAGE1 Matdraw_t Matdraw_lag AGE_t Mat_ETA_AGE Y_tot Matdraw_tot MatAGE_tot...
    Y_t Y_lag meanT stdT K3t Matdraw_eps C A A_tot  C_tot meanA stdA meanC stdC A_t A_lag C_lag ...
    Vect_c_xi  Vect_c_age Vect_c_eta Vect_c_eps Vect_c_a  Vect_a_xi Vect_a_age Vect_a_eta ...
    Vect_a_eps Vect_a_c Vect_a_a Vect_xi_a Vect_xi_age Vect_xi_eta Vect_a1_eta Vect_a1_age...
    Knots_c_xi Knots_c_age Knots_c_eta Knots_c_eps Knots_c_a Knots_a_xi Knots_a_age Knots_a_eta ...
    Knots_a_eps Knots_a_c Knots_a_a Knots_xi_a Knots_xi_age Knots_xi_eta Knots_a1_eta Knots_a1_age...
    AGE_tot M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13 M14 M15 M16 MatXi MatClag Matdraw MatA MatAlag...
    uc Vect_AGE_tot Vect_A_tot Vect_AGE_t Vect_A_lag Vect_Matdraw_lag Vect_Y_lag Vect_C_lag ua ua1 Vect_AGE1 Vect_Matdraw1 unit_interval...
    EDUC YB meanEDUC meanYB stdEDUC stdYB Ybar  meanYbar stdYbar M12a M12b M17   Vect_YB Vect_ED...
    OLSnew OLSnew_e0 OLSnew_eps veps veta veta1 b1 bL b1_e0 bL_e0 b1_eps bL_eps Resqfinal_e0 Resqfinal_eps Resqfinal

%%
% 1. Parameters 
%%%%%%%%%%%%

% parpool('local',20)

peturb_init=0; % Boolean to toggle random perturbation of initial conditions
asset_rule=1; % Boolean to toggle whether asset rule appears in the likelihood / posterior
Nis_scale=50; % Scaling parameter for importance samples
tol_dnom=2; % Resampling tolerance
draws=50; % Length of chains in MH step
sigma=1; % s.d. for the standard normals used in initial conditions
MH_scale=3.5; % Random walk proposal variance as multiple of std(xi) (i.e. proposal variance = MH_scale*sigma_xi)
maxiter=maxitercons; % Total StEM iterations

% Order of hermite polynomials (others will be loaded from income results)
M1=2;
M2=2;
M3=1;
M4=1;
M5=1;
K=2;
M6=2;
M7=1; % This is lagged INCOME (i.e. NOT lagged eta)
M8=0; % This is lagged epsilon - (ALWAYS SET TO 0 IN BASELINE)
M9=1;
M10=1;
M11=0;
M12=1;
M12a=1;
M12b=1;
M13=1;
M14=1;
M15=1;
M16=1;
M17=1;

%%
% 2. Load data 
%%%%%%%%%%%%

% Load earnings parameters

load(savepathincome,'T','Ntau', 'Vectau',...
    'Matdraw_tot', 'Matdraw_lag',...
    'K1', 'K2', 'K2t', 'K3', 'K3t', 'K4', 'K5', 'K6', 'Matdraw', 'Matdraw1',...
    'Resqfinal_e0', 'Resqfinal_eps', 'Resqfinal',...
    'b1', 'bL', 'b1_e0', 'bL_e0', 'b1_eps', 'bL_eps',...
    'OLSnew', 'OLSnew_e0', 'OLSnew_eps', 'veps', 'veta', 'veta1')


% load income data parameters (avoids saving them which takes up lots of
% storage)
load_data(data);

% Load the consumption data
load_cons_data(data);

[Resqinit_cons,Resqinit_a,Resqinit_xi,Resqinit_a1,OLSnew_cons,OLSnew_a,OLSnew_xi,OLSnew_a1,vcons,va,vxi,va1,b1_c,bL_c,b1_a,bL_a,b1_a1,bL_a1,b1_xi,bL_xi]=init_cons_pmcmc_simple();

Resqinit= Resqfinal;
Resqinit_e0 = Resqfinal_e0;
Resqinit_eps=Resqfinal_eps;

Resqnew_cons=zeros((M1+1)*(M2+1)*(M3+1)*(M4+1)*(M5+1),Ntau,maxiter);
Resqnew_a=zeros((M6+1)*(M7+1)*(M8+1)*(M9+1)*(M10+1)*(M11+1),Ntau,maxiter);
Resqnew_xi=zeros((M12+1)*(M12a+1)*(M12b+1),Ntau,maxiter);
Resqnew_a1=zeros((M13+1)*(M14+1)*(M15+1)*(M16+1)*(M17+1),Ntau,maxiter);


% Initialize objects needed for StEM
mat_b=zeros(maxiter,8);
Matdraw=zeros(N,T+1);
Matdraw_store=zeros(N,T+1,maxiter);
vstore=[];
xi_draws=zeros(N,1);
ESS_store=zeros(N,T,maxiter);
Acceptrate_store=zeros(N,draws,maxiter);

% initialize with estimated paramters to speed things up
res = load ('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Cons/Results/REVISION_FINAL.mat',...
    'Resqfinal_cons','Resqfinal_a','Resqfinal_a1','Resqfinal_xi',...
    'OLSnew_cons','OLSnew_a','OLSnew_a1','OLSnew_xi',...
    'vcons','va','va1','vxi','b1_c','bL_c','b1_a','bL_a','b1_a1','bL_a1','b1_xi','bL_xi');

Resqinit_cons = res.Resqfinal_cons;
Resqinit_a = res.Resqfinal_a;
Resqinit_a1 = res.Resqfinal_a1;
Resqinit_xi = res.Resqfinal_xi;
OLSnew_cons = res.OLSnew_cons;
OLSnew_a = res.OLSnew_a;
OLSnew_xi = res.OLSnew_xi;
OLSnew_a1 = res.OLSnew_a1;
vcons = res.vcons;
va = res.va;
vxi = res.vxi;
va1 = res.va1;
b1_c = res.b1_c;
bL_c = res.bL_c;
b1_a = res.b1_a;
bL_a = res.bL_a;
b1_a1 = res.b1_a1;
bL_a1 = res.bL_a1;
b1_xi = res.b1_xi;
bL_xi = res.bL_xi;

% If toggled randomly perturb the initial draws
if peturb_init==1
Resqinit_cons=Resqinit_cons+.1*randn(size(Resqinit_cons));
Resqinit_a=Resqinit_a+.1*randn(size(Resqinit_a));
Resqinit_a1=Resqinit_a1+.1*randn(size(Resqinit_a1));
Resqinit_xi=Resqinit_xi+.1*randn(size(Resqinit_xi));
b1_c=b1_c*2*rand(1);
bL_c=bL_c*2*rand(1);
b1_a=b1_a*2*rand(1);
bL_a=bL_a*2*rand(1);
b1_xi=b1_xi*2*rand(1);
bL_xi=bL_xi*2*rand(1);
end

%%
% 3. Stochastic EM
%%%%%%%%%%%%

for iter=1:maxiter
    disp('Consumption iteration')
    disp(iter)
    disp(datetime)
    %%% PMCMC %%%
    
    Matdraw_final=zeros(N,T+1); % Stores the samples, last column is \xi
    acceptrate=zeros(N,draws); % Stores acceptances at each MH iteration here
    lik_iter=zeros(N,1); % Stores marginal likelihood
    ESS_iter=zeros(N,T); % Stores effective samle sizes
    
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
                end
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
                    % (NB: in baseline asset rule is irrelevant and will cancel from likelihood)
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
    
    
    % Useful high-level summaries of ESS
    ESS_disp=zeros(1,T);
    ESS_disp2=zeros(1,T);
    for tt=1:T
        ESS_temp=ESS_iter(:,tt);
        sample=find(ESS_temp);
        ESS_disp(tt)=mean( ESS_temp(sample));
        ESS_disp2(tt)=nanmean(ESS_temp(sample));
    end
    ESS_disp;
    ESS_disp2;
    ESS_store(:,:,iter)=ESS_iter;
    
    Matdraw=Matdraw_final;
    Acceptrate_store(:,:,iter)=(acceptrate);
    mean(acceptrate);
    
    options.Display ='off';
    warning off
    
    %%% Normalize xi
    
    sample=find(D(:));
    Matdraw_tot=Matdraw(:,1:T);
    Matdraw_tot=Matdraw_tot(:);
    Matdraw_tot=Matdraw_tot(sample);
    
    Matzeta_tot=Matdraw(:,T+1);
    for j=2:T
        Matzeta_tot=[Matzeta_tot;Matdraw(:,T+1)];
    end
    Matzeta_tot=Matzeta_tot(sample);
    
    % Generate regressors - CONSUMPTION function
    Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(Matdraw_tot-meanY)/stdY),uc(:,2),'Uniform',0);
    Vect_Eps_tot=arrayfun(@(x) hermite(x,(Y_tot-Matdraw_tot-meanY)/stdY),uc(:,3),'Uniform',0);
    Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matzeta_tot-meanC)/stdC),uc(:,5),'Uniform',0);
    MatClag = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});
    
    % Estimate via OLS
    b_ols=pinv(MatClag)*C_tot;
    
    % Generate regressors - CONSUMPTION function evaluated at fixed points
    Vect_AGE_temp=arrayfun(@(x) hermite(x,0),uc(:,4),'Uniform',0);
    Vect_A_temp=arrayfun(@(x) hermite(x,0),uc(:,1),'Uniform',0);
    Vect_Matdraw_temp=arrayfun(@(x) hermite(x,0),uc(:,2),'Uniform',0);
    Vect_Eps_temp=arrayfun(@(x) hermite(x,0),uc(:,3),'Uniform',0);
    Vect_Xi_temp=arrayfun(@(x) hermite(x,(Matdraw(:,T+1)-meanC)/stdC),uc(:,5),'Uniform',0);
    MatClag = cat(2,Vect_A_temp{:}).*cat(2,Vect_Matdraw_temp{:}).*cat(2,Vect_Eps_temp{:}).*cat(2,Vect_AGE_temp{:}).*cat(2,Vect_Xi_temp{:});
    
    % Reset xi
    Matdraw(:,T+1) = MatClag*b_ols;
      
    %%% Main M-Step
    
    % Generate regressors - CONSUMPTION function
    sample=find(D(:));
    Matdraw_tot=Matdraw(:,1:T);
    Matdraw_tot=Matdraw_tot(:);
    Matdraw_tot=Matdraw_tot(sample);
    Matzeta_tot=Matdraw(:,T+1);
    for j=2:T
        Matzeta_tot=[Matzeta_tot;Matdraw(:,T+1)];
    end
    Matzeta_tot=Matzeta_tot(sample);
    Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(Matdraw_tot-meanY)/stdY),uc(:,2),'Uniform',0);
    Vect_Eps_tot=arrayfun(@(x) hermite(x,(Y_tot-Matdraw_tot-meanY)/stdY),uc(:,3),'Uniform',0);
    Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matzeta_tot-meanC)/stdC),uc(:,5),'Uniform',0);
    MatClag = cat(2,Vect_A_tot{:}).*cat(2,Vect_Matdraw_tot{:}).*cat(2,Vect_Eps_tot{:}).*cat(2,Vect_AGE_tot{:}).*cat(2,Vect_Xi_tot{:});
    
    % Generate regressors - ASSETS function
    sample=find(D_t);
    i_mean_C_t=Matdraw(:,T+1);
    for tt=3:T
        i_mean_C_t=[i_mean_C_t Matdraw(:,T+1)];
    end
    i_mean_C_t=i_mean_C_t(sample);
    Matdraw_lag=Matdraw(:,1:T-1);
    Matdraw_lag=Matdraw_lag(sample);
    Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(Y_lag-meanY)/stdY),ua(:,2),'Uniform',0);
    Vect_Xi_lag=arrayfun(@(x) hermite(x,(i_mean_C_t-meanC)/stdC),ua(:,6),'Uniform',0);
    Vect_Y_lag=arrayfun(@(x) hermite(x,(Y_lag-Matdraw_lag-meanY)/stdY),ua(:,3),'Uniform',0);
    MatAlag = cat(2,Vect_A_lag{:}).*cat(2,Vect_C_lag{:}).*cat(2,Vect_Y_lag{:}).*cat(2,Vect_Matdraw_lag{:}).*cat(2,Vect_AGE_t{:}).*cat(2,Vect_Xi_lag{:});
    
    % Generate regressors - ASSETS function
    Matdraw1=zeros(N,1);
    for iii=1:N
        Matdraw1(iii)=Matdraw(iii,ENT(iii));
    end
    Vect_Matdraw1=arrayfun(@(x) hermite(x,(Matdraw1-meanY)/stdY),ua1(:,1),'Uniform',0);
    Vect_Xi1=arrayfun(@(x) hermite(x,(Matdraw(:,T+1)-meanC)/stdC),ua1(:,3),'Uniform',0);
    MatA =cat(2,Vect_Matdraw1{:}).*cat(2,Vect_AGE1{:}).*cat(2,Vect_Xi1{:}).*cat(2,Vect_YB{:}).*cat(2,Vect_ED{:});

    % Quantile regressions
    for jtau=1:Ntau
        tau=Vectau(jtau);
        Resqnew_cons(:,jtau,iter)=fminunc(@wqregk_cons_age,Resqinit_cons(:,jtau),options);
        Resqnew_xi(:,jtau,iter)=fminunc(@wqregk_xi_age,Resqinit_xi(:,jtau),options);
        
        % Only update asset rule if toggled
        if asset_rule==1
        Resqnew_a(:,jtau,iter)=fminunc(@wqregk_a_age,Resqinit_a(:,jtau),options);
        Resqnew_a1(:,jtau,iter)=fminunc(@wqregk_a1_age,Resqinit_a1(:,jtau),options);
        end
        
    end
    
    % Update MH variance
    beta=pinv(MatXi)*Matdraw(:,T+1);
    vxi=var(Matdraw(:,T+1)-MatXi*beta);
    
    % Taste heterogeneity: Laplace parameters
    Vect1=C_tot-MatClag*Resqnew_cons(:,1,iter);
    Vect2=C_tot-MatClag*Resqnew_cons(:,Ntau,iter);
    b1_c=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_c=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    Resqinit_cons=Resqnew_cons(:,:,iter);

    Vect1=Matdraw(:,T+1)-MatXi*Resqnew_xi(:,1,iter);
    Vect2=Matdraw(:,T+1)-MatXi*Resqnew_xi(:,Ntau,iter);
    b1_xi=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_xi=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    Resqinit_xi=Resqnew_xi(:,:,iter);
    
    % If toggled estimate asset rule parameters
    if asset_rule==1
    Vect1=A_t-MatAlag*Resqnew_a(:,1,iter);
    Vect2=A_t-MatAlag*Resqnew_a(:,Ntau,iter);
    b1_a=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_a=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    Vect1=A1-MatA*Resqnew_a1(:,1,iter);
    Vect2=A1-MatA*Resqnew_a1(:,Ntau,iter);
    b1_a1=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_a1=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    Resqinit_a=Resqnew_a(:,:,iter);
    Resqinit_a1=Resqnew_a1(:,:,iter);
    end
    
    % Store other generated values
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
    Matdraw_store(:,:,iter) = Matdraw;
    
end

% Prepare results for saving

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

disp('Saving consumption results')

b1_c=mean(mat_b((maxiter/2):maxiter,1))
bL_c=mean(mat_b((maxiter/2):maxiter,2))
b1_xi=mean(mat_b((maxiter/2):maxiter,3))
bL_xi=mean(mat_b((maxiter/2):maxiter,4))
b1_a=mean(mat_b((maxiter/2):maxiter,5))
bL_a=mean(mat_b((maxiter/2):maxiter,6))
b1_a1=mean(mat_b((maxiter/2):maxiter,7))
bL_a1=mean(mat_b((maxiter/2):maxiter,8))


save(savepathconsumption,'data','b1_c','bL_c','b1_xi','bL_xi','b1_a','bL_a','b1_a1','bL_a1','Resqfinal_cons','Resqfinal_xi','Resqfinal_a','Resqfinal_a1','mat_lik','Matdraw_store')
end
