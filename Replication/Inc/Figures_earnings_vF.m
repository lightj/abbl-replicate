% ABBL (2021)
% Plot main earnings results

cd '/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/Inc'

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

load ('/home/jdlight/ABBL - PMCMC/JOE_codes/Results/211015_smc_n.mat')

figs(1)=figure;
plot(mat_lik)
xlim([0, maxiter])
ylabel('Likelihood','FontSize',9)
xlabel('Iteration','FontSize',9)

Resqtrue=Resqfinal;
Resqtrue_e0=Resqfinal_e0;
Resqtrue_eps=Resqfinal_eps;

b1true=b1;
bLtrue=bL;
b1true_e0=b1_e0;
bLtrue_e0=bL_e0;
b1true_eps=b1_eps;
bLtrue_eps=bL_eps;

Resqinit=Resqfinal;
Resqinit_e0=Resqfinal_e0;
Resqinit_eps=Resqfinal_eps;

%%
% 2. Simulate posterior draws
%%%%%%%%%%%%

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

Mateta_true=Matdraw_final;
Mateps_true=Y-Mateta_true;
Ytilde=Mateta_true+Mateps_true;

Y1=zeros(N,1);
for iii=1:N
    Y1(iii)=Y(iii,ENT(iii));
end

%%
% 3. Figures in paper
%%%%%%%%%%%%

Ks=3;
Mat1=[];
for kk1=0:Ks
    Mat1=[Mat1 hermite(kk1,(Y_lag-meanY)/stdY)];
end

ResP_data=zeros((Ks+1),Ntau);
for jtau=1:Ntau
    tau=Vectau(jtau);
    ResP_data(:,jtau)=rq(Mat1,Y_t(:),tau);
end
Mat2=zeros(size(Y_lag,1),1);
for kk1=1:Ks
    Mat2=[Mat2 kk1*hermite(kk1-1,(Y_lag-meanY)/stdY)./stdY];
end
mean(Mat2*ResP_data)

Vect=quantile(Y_lag,Vectau);
Mat3=zeros(Ntau,1);
for kk1=1:Ks    
    Mat3=[Mat3 kk1*hermite(kk1-1,(Vect-meanY)/stdY)./stdY];
end

Mat3*ResP_data

Pers_y = Mat3*ResP_data;
figs(1)=figure;
set(figs(1), 'Position', [10 10 400 400]);
surf(Vectau,Vectau,Pers_y)
xlabel('Income shock','FontSize',11)
ylabel('Initial income','FontSize',11)
zlabel('Persistence','FontSize',11)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1.2])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1.2))
xtickangle(0)
ytickangle(0)
text(0.2,0.9,0.35,['\mu = ' num2str(mean(Pers_y(:)),2)])
text(0.2,0.9,0.25,['\sigma = ' num2str(std(Pers_y(:)),2)])
caxis([0 1])
print(figs(1),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/figs/abb_data_%d.eps',model),'-depsc');

ages=unique(AGE_t);
prob_age=zeros(size(ages,1),1);
for aaa=1:size(ages,1)
    prob_age(aaa)=mean(AGE_t(:)==ages(aaa));
end

sample=find(D_t(:));
Vect=Mateta_true(:,1:T-1);
Vect=Vect(:);
Vect=Vect(sample);

Vect=quantile(Vect(:),Vectau);

Pers_eta=zeros(Ntau,Ntau);

for aaa=1:size(ages,1)
    
    Mat=zeros(Ntau,(K2+1)*(K2t+1));
    for kk1=1:K1
        for kk2=0:K2
            for kk2t=0:K2t
                Mat=[Mat kk1*hermite(kk1-1,(Vect(:)-meanY)/stdY)./stdY.*...
                    hermite(kk2,(ages(aaa)-meanAGE)/stdAGE)...
                    .*hermite(kk2t,4)];
            end
        end
    end
    Pers_eta=Pers_eta+prob_age(aaa)*Mat*Resqtrue;
end

figs(2)=figure;
set(figs(2), 'Position', [10 10 400 400]);
surf(Vectau,Vectau,Pers_eta)
xlabel('Income shock','FontSize',11)
ylabel('Initial income','FontSize',11)
zlabel('Persistence','FontSize',11)
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'zlim',[0 1.2])
set(gca,'xtick',(0:0.2:1))
set(gca,'ytick',(0:0.2:1))
set(gca,'ztick',(0:0.2:1.2))
xtickangle(0)
ytickangle(0)
text(0.2,0.9,0.35,['\mu = ' num2str(mean(Pers_eta(:)),2)])
text(0.2,0.9,0.25,['\sigma = ' num2str(std(Pers_eta(:)),2)])
caxis([0 1])
print(figs(2),sprintf('/home/jdlight/ABBL - PMCMC/JOE_codes/SMC/figs/abb_pers_%d.eps',model),'-depsc');

