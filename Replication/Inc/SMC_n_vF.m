% StEM estimation using SMC methods to estimate the model in Arellano, Blundell and Bonhomme (2017)
% This version: September 2021
% Questions: jdlight@uchicago.edu
  
clear all
clc;

cd '/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/Inc'

global T N K1 K2 K3 K4 K5 Ntau Vectau tau ...
       D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
       Matdraw1 MatAGE1 Matdraw_t Matdraw_lag AGE_t Mat_ETA_AGE Y_tot Matdraw_tot MatAGE_tot...
       Y_t Y_lag meanT stdT K3t Matdraw_eps Time_t K2t K6 EDUC YB meanEDUC meanYB stdEDUC stdYB

% parpool('local',56)

data=load ('/home/jdlight/Data/data4est_vFINAL_WITH_DUAL_RESTRICTION.out');

%%% Set parameters %%%

T=7; % Time periods
Nis=5000; % No. of importance samples in E-step
tol=Nis/2; % Resampling tolerance
maxiter=500; % # Stochastic EM iterations
rng('shuffle') % Generate seed
Ntau=11; % Set grid
Vectau=(1/(Ntau+1):1/(Ntau+1):Ntau/(Ntau+1))';

K1=3; % Hermite eta on lag eta
K2=2; % Hermite eta on age
K2t=0;
K3=2; % Hermite eps on age
K3t=0; % Hermite eps on calendar year
K4=2; % Hermite eta1 on age
K5=1; % Hermite eta1 on EDUC
K6=1; % Hermite eta1 on YB

% Load the data
load_data(data);

% Generate initial values
[Resqinit,Resqinit_e0,Resqinit_eps,OLSnew,OLSnew_e0,OLSnew_eps,veta1,veta,veps,b1,bL,b1_e0,bL_e0,b1_eps,bL_eps]=init();

Resqnew=zeros((K1+1)*(K2+1)*(K2t+1),Ntau,maxiter);
Resqnew_e0=zeros((K4+1)*(K5+1)*(K6+1),Ntau,maxiter);
Resqnew_eps=zeros((K3+1)*(K3t+1),Ntau,maxiter);

% b1=10;
% bL=10;
% b1_e0=10;
% bL_e0=10;
% b1_eps=10;
% bL_eps=10;

b1init=b1;
bLinit=bL;
b1init_e0=b1_e0;
bLinit_e0=bL_e0;
b1init_eps=b1_eps;
bLinit_eps=bL_eps;

Resqinit_eps=Resqinit_eps-mean(Resqinit_eps')'*ones(1,Ntau);
Resqinit_eps(1,:)=Resqinit_eps(1,:)-((1-Vectau(Ntau))/bL_eps-Vectau(1)/b1_eps)*ones(1,Ntau);
    
% Storage matrices used in the routine
mat_b=zeros(maxiter,6);
mat_lik=zeros(maxiter,1);
ESS_store=zeros(N,T,maxiter);

OLSstore_eta=[];
OLSstore_eta1=[];
OLSstore_eps=[];
var_store=[];
Matdraw=zeros(N,T);

for iter=1:maxiter
    iter
    
    %E step
    
    %%%%%%%%%%%%%%% Adaptive SMC %%%%%%%%%%%%%%%
    
    Matdraw_final=zeros(N,T);
    lik_iter=zeros(N,1);
    ESS_iter=zeros(N,T);
    parfor iii=1:N
        
        ESS_i=zeros(1,T);
        T1=ENT(iii);
        TT=EXT(iii);
        
        
        i_Y=ones(Nis,1)*Y(iii,:);
        i_AGE=ones(Nis,1)*AGE(iii,:);
        i_time=ones(Nis,1)*(1:1:T);
        
        i_MatAGE_tot=[];
        for kk3=0:K3
            for kk3t=0:K3t
            i_MatAGE_tot=[i_MatAGE_tot hermite(kk3,(i_AGE(:)-meanAGE)/stdAGE).*hermite(kk3t,(i_time(:)-meanT)/stdT)];
            end
        end
        
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
        
        % Define useful matrices
        particle_draw=zeros(Nis,T);
        mat_resamp=zeros(T,1);
        
        
        for ttt=T1:TT
            if ttt==T1
                alpha=veta1/(veta1 + veps);
                Mvar=2/((1/veta1)+(1/veps));
                
                Mmu=(1-alpha)*(i_MatAGE1(1,:)*OLSnew_e0)...
                    + alpha*i_Y(1,ttt);
                
                % Generate the importance samples and their weights
                particle_draw(:,ttt) = Mmu + mvnrnd(0,Mvar,Nis);
                weights=IS_weight(i_Y,particle_draw,i_MatAGE_tot,i_MatAGE1,i_AGE,ttt,...
                    Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,...
                    b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,K1,K2,meanAGE,...
                    stdAGE,meanY,stdY,Nis,T1,[]);
                
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
                            .*hermite(kk2,(i_AGE(:,ttt)-meanAGE)/stdAGE).*hermite(kk2t,ttt)];
                        end
                    end
                end
                
               Mmu=(1-alpha)*(i_Matdraw_lag*OLSnew)...
                    + alpha.*i_Y(:,ttt);
                
                
                % Generate the importance samples and their weights
                particle_draw(:,ttt) = Mmu + mvnrnd(0,Mvar,Nis);
                weights=weights...
                    .*IS_weight(i_Y,particle_draw,i_MatAGE_tot,i_MatAGE1,i_AGE,ttt,...
                    Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,...
                    b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,K1,K2,meanAGE,...
                    stdAGE,meanY,stdY,Nis,T1,i_Matdraw_lag);
                
                weights=weights./mvnpdf(particle_draw(:,ttt)-Mmu,0,Mvar);
                % weights(isnan(weights))=0;
                lik=lik*mean(weights);
                
            end
            lik_iter(iii)=lik;
            % Self-normalization
            weights=weights./sum(weights);
            
            % Effective sample size
            ESS= 1/sum(weights.^2);
            ESS_i(ttt)=ESS;
%           
            % Adaptive resampling (tolerance set by tol)
            if ttt<TT && ESS<tol
                particle_draw = datasample(particle_draw,Nis,'weights',weights);
                weights=ones(Nis,1)*(1/Nis);
            elseif ttt==TT
                Matdraw_old=datasample(particle_draw,1,'weights',weights);
            end
            newObj=lik;
            
        end
        
        
        Matdraw_final(iii,:) = Matdraw_old;
        
        ESS_iter(iii,:)=ESS_i;
        
    end
    
    ESS_store(:,:,iter)=ESS_iter;
    ESS_mean=zeros(T,1);
    for tt=1:T
    ESS_mean(tt)=mean(nonzeros(ESS_iter(:,tt)));    
    end
    ESS_mean'    %Last draws of the chain will be the fixed associated with our data.
    Matdraw=Matdraw_final;

    options.Display ='off';
    warning off
    
    sample=find(D_t);
    
    Matdraw_t=Matdraw(:,2:T);
    Matdraw_t=Matdraw_t(:);
    Matdraw_t=Matdraw_t(sample);
    
    Matdraw_lag=Matdraw(:,1:T-1);
    Matdraw_lag=Matdraw_lag(:);
    Matdraw_lag=Matdraw_lag(sample);
    
    Mat_ETA_AGE=[];
    for kk1=0:K1
        for kk2=0:K2
            for kk2t=0:K2t
            Mat_ETA_AGE=[Mat_ETA_AGE hermite(kk1,(Matdraw_lag-meanY)/stdY).*hermite(kk2,(AGE_t-meanAGE)/stdAGE).*hermite(kk2t,Time_t)];
            end
        end
    end
    
    Matdraw1=zeros(N,1);
    for ii=1:N
       Matdraw1(ii)=Matdraw(ii,ENT(ii)); 
    end
    
    sample=find(D(:));
    Matdraw_tot=Matdraw(:);
    Matdraw_tot=Matdraw_tot(sample);
    
    OLSnew_e0=pinv(MatAGE1)*Matdraw1;
    OLSnew=pinv(Mat_ETA_AGE)*Matdraw_t;
    OLSnew_eps=pinv(MatAGE_tot)*(Y_tot-Matdraw_tot);
    
    Matdraw_eps=Y_tot-Matdraw_tot;
    veta1=var(Matdraw1-MatAGE1*OLSnew_e0);
    veta=var(Matdraw_t-Mat_ETA_AGE*OLSnew);
    veps=var(Matdraw_eps);
    
    
    for jtau=1:Ntau
        
        tau=Vectau(jtau);
        
        Resqnew_e0(:,jtau,iter)=fminunc(@wqregk_e0_age,Resqinit_e0(:,jtau),options);
        Resqnew(:,jtau,iter)=fminunc(@wqregk_pt_age,Resqinit(:,jtau),options);
        Resqnew_eps(:,jtau,iter)=fminunc(@wqregk_eps_age,Resqinit_eps(:,jtau),options);
        
    end
        
    
    %% Normalization
     Resqnew_eps(:,:,iter)=Resqnew_eps(:,:,iter)-mean(Resqnew_eps(:,:,iter)')'*ones(1,Ntau);
     Resqnew_eps(1,:,iter)=Resqnew_eps(1,:,iter)-((1-Vectau(Ntau))/bL_eps-Vectau(1)/b1_eps)*ones(1,Ntau);
    

    warning on
    
    % Laplace parameters: draws
    Vect1=Matdraw_eps-MatAGE_tot*Resqnew_eps(:,1,iter);
    Vect2=Matdraw_eps-MatAGE_tot*Resqnew_eps(:,Ntau,iter);
    b1_eps=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_eps=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    
    
    Vect1=Matdraw1-MatAGE1*Resqnew_e0(:,1,iter);
    Vect2=Matdraw1-MatAGE1*Resqnew_e0(:,Ntau,iter);
    b1_e0=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL_e0=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    
    Vect1=Matdraw_t-Mat_ETA_AGE*Resqnew(:,1,iter);
    Vect2=Matdraw_t-Mat_ETA_AGE*Resqnew(:,Ntau,iter);
    b1=-sum(Vect1<=0)/sum(Vect1.*(Vect1<=0));
    bL=sum(Vect2>=0)/sum(Vect2.*(Vect2>=0));
    
    
    
    
    
    % Criterion
    
    
    Resqinit_e0=Resqnew_e0(:,:,iter)
    Resqinit_eps=Resqnew_eps(:,:,iter)
    Resqinit=Resqnew(:,:,iter)
    
    OLSstore_eta=[OLSstore_eta OLSnew];
    OLSstore_eta1=[OLSstore_eta1 OLSnew_e0];
    OLSstore_eps=[OLSstore_eps OLSnew_eps];
    var_store=[var_store [veta ; veta1;veps]];
    mat_b(iter,1)=b1;
    mat_b(iter,2)=bL;
    mat_b(iter,3)=b1_e0;
    mat_b(iter,4)=bL_e0;
    mat_b(iter,5)=b1_eps;
    mat_b(iter,6)=bL_eps;
    
    
    
    % Quick computation of persistence
    sample=find(D_t(:));
    Vect=Matdraw(:,1:T-1);
    Vect=Vect(:);
    Vect=Vect(sample);
    
    Mat=[];
    for kk1=0:K1
        for kk2=0:K2
            for kk2t=0:K2t
            if kk1<1
                Mat=[Mat zeros(size(Vect,1),1)];
            else
                Mat=[Mat kk1*hermite(kk1-1,(Vect(:)-meanY)/stdY)./stdY.*hermite(kk2,(AGE_t-meanAGE)/stdAGE).*hermite(kk2t,4)];
            end
            end
        end
    end
    
    mean(Mat*Resqinit)
    
    Vect=quantile(Vect(:),Vectau);
    age_ref=meanAGE;
    Mat=zeros(Ntau,(K2+1)*(K2t+1));
    for kk1=1:K1
        for kk2=0:K2
            for kk2t=0:K2t
            Mat=[Mat kk1*hermite(kk1-1,(Vect(:)-meanY)/stdY)./stdY.*hermite(kk2,(age_ref-meanAGE)/stdAGE).*hermite(kk2t,4)];
            end
        end
    end
    
    
    Mat*Resqinit
    
    % Likelihood
    mat_lik(iter)=mean(log(lik_iter));
    mat_lik(iter)   
    
end

Resqfinal=zeros((K1+1)*(K2+1)*(K2t+1),Ntau);
for jtau=1:Ntau
    for p=1:(K1+1)*(K2+1)*(K2t+1)
        Resqfinal(p,jtau)=mean(Resqnew(p,jtau,(maxiter/2):maxiter));
    end
end

Resqfinal_e0=zeros((K4+1)*(K5+1)*(K6+1),Ntau);
for jtau=1:Ntau
    for p=1:(K4+1)*(K5+1)*(K6+1)
        Resqfinal_e0(p,jtau)=mean(Resqnew_e0(p,jtau,(maxiter/2):maxiter));
    end
end


Resqfinal_eps=zeros((K3+1)*(K3t+1),Ntau);
for jtau=1:Ntau
    for p=1:(K3+1)*(K3t+1)
        Resqfinal_eps(p,jtau)=mean(Resqnew_eps(p,jtau,(maxiter/2):maxiter));
    end
end

Resqinit=Resqfinal;
Resqinit_e0=Resqfinal_e0;
Resqinit_eps=Resqfinal_eps;

b1=mean(mat_b((maxiter/2):maxiter,1))
bL=mean(mat_b((maxiter/2):maxiter,2))
b1_e0=mean(mat_b((maxiter/2):maxiter,3))
bL_e0=mean(mat_b((maxiter/2):maxiter,4))
b1_eps=mean(mat_b((maxiter/2):maxiter,5))
bL_eps=mean(mat_b((maxiter/2):maxiter,6))

OLSnew=mean(OLSstore_eta(:,(maxiter/2):maxiter),2);
OLSnew_e0=mean(OLSstore_eta1(:,(maxiter/2):maxiter),2);
OLSnew_eps=mean(OLSstore_eps(:,(maxiter/2):maxiter),2);
veta=mean(var_store(1,(maxiter/2):maxiter),2);
veta1=mean(var_store(2,(maxiter/2):maxiter),2);
veps=mean(var_store(3,(maxiter/2):maxiter),2);


save '/home/jdlight/ABBL - PMCMC/JOE_codes/Results/211015_smc_n.mat'

%% Visualize persistence in chain


pers_chain=zeros(maxiter,1);

for ix =1:maxiter
ix
    
Resqfinal=Resqnew(:,:,ix);


Mat=[];
for kk1=0:K1
    for kk2=0:K2
        for kk2t=0:K2t
            if kk1<1
                Mat=[Mat zeros(size(AGE_t,1),1)];
            else
                Mat=[Mat kk1*hermite(kk1-1,(Matdraw_lag-meanY)/stdY)./stdY...
                    .*hermite(kk2,(AGE_t-meanAGE)/stdAGE)];
            end
        end
    end
end
pers_chain(ix)=mean(mean(Mat*Resqfinal));



end

figs(1)=figure;
set(figs(1), 'Position', [10 10 800 300]);
plot(pers_chain)
hold on
ylim([0.75 1])
xlim([0 2000])
xlabel('Iteration')
ylabel('Average persistence')

figs(2)=figure;
set(figs(2), 'Position', [10 10 800 300]);
autocorr(pers_chain(1:maxiter))
