function boot_data = simulate_panel()

% ABBL (2021)
% Code to generate the simulated panels used in the parametric bootstrap

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
model = 1;

%%
% 1. Load parameters and data
%%%%%%%%%%%%

load ('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Cons/Results/REVISION_FINAL.mat')

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

%%% UPPER AND LOWER BOUNDS USED WHEN SIMULATING

Yub=10*max(max(Y));
temp=Y;
temp(temp==0)=NaN;
Ylb=10*min(nanmin(temp));

Aub=10*max(max(A));
temp=A;
temp(temp==0)=NaN;
Alb=10*min(nanmin(temp));

Cub=10*max(max(C));
temp=C;
temp(temp==0)=NaN;
Clb=10*min(nanmin(temp));

%%% PLACEHOLDERS FOR THE SIMULATIONS

Mateta_true=zeros(N,T);
Mateps_true=zeros(N,T);
Matc_true=zeros(N,T);
Mata_true=zeros(N,T);

Mateta_1=zeros(N,1);
Matxi_true = zeros(N,1);
Ytilde=zeros(N,T);

%%
% 2. Simulate panel
%%%%%%%%%%%%

for iii=1:N
    
    Mateta_true_i=zeros(1,T);
    Mateps_true_i=zeros(1,T);
    Mata_true_i = zeros(1,T);
    Matc_true_i = zeros(1,T);
    
    truncated=1;
    while truncated==1
        
        T1=ENT(iii);
        TT=EXT(iii);
        
        % Proposal, eta_1
        V_draw=unifrnd(0,1,1,1);
        
        %First quantile
        Mateta_true_i(T1)=(MatAGE1(iii,:)*Resqtrue_e0(:,1)).*(V_draw<=Vectau(1));
        for jtau=2:Ntau
            Mateta_true_i(T1)=Mateta_true_i(T1)+((MatAGE1(iii,:)*(Resqtrue_e0(:,jtau)-Resqtrue_e0(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
                (V_draw-Vectau(jtau-1))+MatAGE1(iii,:)*Resqtrue_e0(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
        end
        %Last quantile.
        Mateta_true_i(T1)=Mateta_true_i(T1)+(MatAGE1(iii,:)*Resqtrue_e0(:,Ntau)).*(V_draw>Vectau(Ntau));
        
        % Laplace tails
        Mateta_true_i(T1)=Mateta_true_i(T1)+((1/(b1true_e0)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
            -(1/bLtrue_e0*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
        
        
        % Proposal, eta_t
        for tt=T1:TT-1
            
            Mat=zeros(1,(K1+1)*(K2+1));
            for kk1=0:K1
                for kk2=0:K2
                    Mat(:,kk1*(K2+1)+kk2+1)=hermite(kk1,(Mateta_true_i(tt)-meanY)/stdY).*hermite(kk2,(AGE(iii,tt+1)-meanAGE)/stdAGE);
                end
            end
            
            V_draw=unifrnd(0,1,1,1);
            
            %First quantile
            Mateta_true_i(tt+1)=(Mat*Resqtrue(:,1)).*(V_draw<=Vectau(1));
            for jtau=2:Ntau
                Mateta_true_i(tt+1)=Mateta_true_i(tt+1)+...
                    ((Mat*Resqtrue(:,jtau)-Mat*Resqtrue(:,jtau-1))/...
                    (Vectau(jtau)-Vectau(jtau-1)).*...
                    (V_draw-Vectau(jtau-1))+Mat*Resqtrue(:,jtau-1)).*...
                    (V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
            end
            %Last quantile.
            Mateta_true_i(tt+1)=Mateta_true_i(tt+1)+(Mat*Resqtrue(:,Ntau)).*...
                (V_draw>Vectau(Ntau));
            
            % Laplace tails
            Mateta_true_i(tt+1)=Mateta_true_i(tt+1)+((1/(b1true)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
                -(1/bLtrue*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
            
            
        end
        Mateta_1(iii)=Mateta_true_i(T1);
        
        Mateta_true(iii,:)=Mateta_true_i;
        
        for tt=T1:TT
            % Proposal, eps
            V_draw=unifrnd(0,1,1,1);
            
            i_MatAGE_tot=[];
            for kk3=0:K3
                for kk3t=0:K3t
                    i_MatAGE_tot=[i_MatAGE_tot hermite(kk3,(AGE(iii,tt)'-meanAGE)/stdAGE).*hermite(kk3t,(tt-meanT)/stdT)];
                end
            end
            
            %First quantile
            Mateps_true_i(tt)=(i_MatAGE_tot*Resqtrue_eps(:,1)).*(V_draw<=Vectau(1));
            for jtau=2:Ntau
                Mateps_true_i(tt)=Mateps_true_i(tt)+((i_MatAGE_tot*(Resqtrue_eps(:,jtau)-Resqtrue_eps(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
                    (V_draw-Vectau(jtau-1))+i_MatAGE_tot*Resqtrue_eps(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
            end
            %Last quantile.
            Mateps_true_i(tt)=Mateps_true_i(tt)+(i_MatAGE_tot*Resqtrue_eps(:,Ntau)).*(V_draw>Vectau(Ntau));
            
            % Laplace tails
            Mateps_true_i(tt)=Mateps_true_i(tt)+((1/(b1true_eps)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
                -(1/bLtrue_eps*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
            
        end
        
        Mateps_true(iii,:)=Mateps_true_i;
        
        
        Ytilde_i= Mateta_true_i + Mateps_true_i;
        Ytilde_bar_i = mean(Ytilde_i(T1:TT));
        
        i_MatXi=[];
        for kk12=0:M12
            for kk13=0:M12a
                for kk14=0:M12b
                    i_MatXi=[i_MatXi hermite(kk12,(EDUC(iii)-meanEDUC)/stdEDUC)...
                        .*hermite(kk13,(YB(iii)-meanYB)/stdYB)...
                        .*hermite(kk14,(Ytilde_bar_i-meanYbar)/stdYbar)];
                end
            end
        end
        
        % Proposal, xi
        V_draw=unifrnd(0,1,1,1);
        
        %First quantile
        Matxi_true_i=(i_MatXi*Resqtrue_xi(:,1)).*(V_draw<=Vectau(1));
        
        for jtau=2:Ntau
            Matxi_true_i=Matxi_true_i+((i_MatXi*(Resqtrue_xi(:,jtau)-Resqtrue_xi(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
                (V_draw-Vectau(jtau-1))+i_MatXi*Resqtrue_xi(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
        end
        
        %Last quantile.
        Matxi_true_i=Matxi_true_i+(i_MatXi*Resqtrue_xi(:,Ntau)).*(V_draw>Vectau(Ntau));
        
        % Laplace tails
        Matxi_true_i=Matxi_true_i+((1/(b1true_xi)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
            -(1/bLtrue_xi*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
        
        for ttt=T1:TT
            
            V_draw=unifrnd(0,1,1,1);
            
            if ttt==T1
                
                i_Vect_AGE1=arrayfun(@(x) hermite(x,(AGE(iii,ttt)-meanAGE)/stdAGE),ua1(:,2),'Uniform',0);
                i_Vect_Matdraw1=arrayfun(@(x) hermite(x,(Mateta_true_i(ttt)-meanY)/stdY),ua1(:,1),'Uniform',0);
                i_Vect_Xi1=arrayfun(@(x) hermite(x,(Matxi_true_i-meanC)/stdC),ua1(:,3),'Uniform',0);
                i_Vect_YB=arrayfun(@(x) hermite(x,(YB(iii)-meanYB)/stdYB),ua1(:,4),'Uniform',0);
                i_Vect_ED=arrayfun(@(x) hermite(x,(EDUC(iii)-meanEDUC)/stdEDUC),ua1(:,5),'Uniform',0);
                
                i_MatA =cat(2,i_Vect_Matdraw1{:}).*cat(2,i_Vect_AGE1{:}).*cat(2,i_Vect_Xi1{:}).*cat(2,i_Vect_YB{:}).*cat(2,i_Vect_ED{:});
                
                %First quantile
                Mata_true_i(ttt)=(i_MatA*Resqtrue_a1(:,1)).*(V_draw<=Vectau(1));
                for jtau=2:Ntau
                    Mata_true_i(ttt)=Mata_true_i(ttt)+((i_MatA*(Resqtrue_a1(:,jtau)-Resqtrue_a1(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
                        (V_draw-Vectau(jtau-1))+i_MatA*Resqtrue_a1(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
                end
                %Last quantile.
                Mata_true_i(ttt)=Mata_true_i(ttt)+(i_MatA*Resqtrue_a1(:,Ntau)).*(V_draw>Vectau(Ntau));
                
                % Laplace tails
                Mata_true_i(tt)=Mata_true_i(tt)+((1/(b1true_a1)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
                    -(1/bLtrue_a1*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
                
            elseif ttt>T1
                
                i_Vect_AGE_t=arrayfun(@(x) hermite(x,(AGE(iii,ttt)-meanAGE)/stdAGE),ua(:,4),'Uniform',0);
                i_Vect_A_lag=arrayfun(@(x) hermite(x,(Mata_true_i(ttt-1)-meanA)/stdA),ua(:,1),'Uniform',0);
                i_Vect_Matdraw_lag=arrayfun(@(x) hermite(x,(Ytilde_i(ttt-1)-meanY)/stdY),ua(:,2),'Uniform',0);
                i_Vect_Y_lag=arrayfun(@(x) hermite(x,(Ytilde_i(ttt-1)-Mateta_true_i(ttt-1)-meanY)/stdY),ua(:,3),'Uniform',0);
                i_Vect_C_lag=arrayfun(@(x) hermite(x,(Matc_true_i(ttt-1)-meanC)/stdC),ua(:,5),'Uniform',0);
                i_Vect_Xi_lag=arrayfun(@(x) hermite(x,(Matxi_true_i-meanC)/stdC),ua(:,6),'Uniform',0);
                
                i_MatA = cat(2,i_Vect_A_lag{:}).*cat(2,i_Vect_C_lag{:}).*cat(2,i_Vect_Y_lag{:}).*cat(2,i_Vect_Matdraw_lag{:}).*cat(2,i_Vect_AGE_t{:}).*cat(2,i_Vect_Xi_lag{:});
                
                %First quantile
                Mata_true_i(ttt)=(i_MatA*Resqtrue_a(:,1)).*(V_draw<=Vectau(1));
                for jtau=2:Ntau
                    Mata_true_i(ttt)=Mata_true_i(ttt)+((i_MatA*(Resqtrue_a(:,jtau)-Resqtrue_a(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
                        (V_draw-Vectau(jtau-1))+i_MatA*Resqtrue_a(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
                end
                %Last quantile.
                Mata_true_i(ttt)=Mata_true_i(ttt)+(i_MatA*Resqtrue_a(:,Ntau)).*(V_draw>Vectau(Ntau));
                
                % Laplace tails
                Mata_true_i(ttt)=Mata_true_i(ttt)+((1/(b1true_a)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
                    -(1/bLtrue_a*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
                
            end
            i_Vect_AGE_tot=arrayfun(@(x) hermite(x,(AGE(iii,ttt)-meanAGE)/stdAGE),uc(:,4),'Uniform',0);
            i_Vect_A_tot=arrayfun(@(x) hermite(x,(Mata_true_i(ttt)-meanA)/stdA),uc(:,1),'Uniform',0);
            i_Vect_Matdraw_tot=arrayfun(@(x) hermite(x,(Mateta_true_i(ttt)-meanY)/stdY),uc(:,2),'Uniform',0);
            i_Vect_Eps_tot=arrayfun(@(x) hermite(x,(Ytilde_i(ttt)-Mateta_true_i(ttt)-meanY)/stdY),uc(:,3),'Uniform',0);
            i_Vect_Xi_tot=arrayfun(@(x) hermite(x,(Matxi_true_i-meanC)/stdC),uc(:,5),'Uniform',0);
            i_MatC = cat(2,i_Vect_A_tot{:}).*cat(2,i_Vect_Matdraw_tot{:}).*cat(2,i_Vect_Eps_tot{:}).*cat(2,i_Vect_AGE_tot{:}).*cat(2,i_Vect_Xi_tot{:});
            
            V_draw=unifrnd(0,1,1,1);
            
            %First quantile
            Matc_true_i(ttt)=(i_MatC*Resqtrue_c(:,1)).*(V_draw<=Vectau(1));
            for jtau=2:Ntau
                Matc_true_i(ttt)=Matc_true_i(ttt)+((i_MatC*(Resqtrue_c(:,jtau)-Resqtrue_c(:,jtau-1)))/(Vectau(jtau)-Vectau(jtau-1)).*...
                    (V_draw-Vectau(jtau-1))+i_MatC*Resqtrue_c(:,jtau-1)).*(V_draw>Vectau(jtau-1)).*(V_draw<=Vectau(jtau));
            end
            %Last quantile.
            Matc_true_i(ttt)=Matc_true_i(ttt)+(i_MatC*Resqtrue_c(:,Ntau)).*(V_draw>Vectau(Ntau));
            
            % Laplace tails
            Matc_true_i(ttt)=Matc_true_i(ttt)+((1/(b1true_c)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
                -(1/bLtrue_c*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
            
            
        end
        Ytest=Mateta_true_i(T1:TT) + Mateps_true_i(T1:TT);
        test_uby = max(Ytest);
        test_lby = min(Ytest);
        
        test_ub = max(Mata_true_i(T1:TT));
        test_lb = min(Mata_true_i(T1:TT));
        test_ub2 = max(Matc_true_i(T1:TT));
        test_lb2 = min(Matc_true_i(T1:TT));
        
        truncated = 1-(test_ub<=Aub)*(test_lb>=Alb)*(test_ub2<=Cub)*(test_lb2>=Clb)*(test_uby<=Yub)*(test_lby>=Ylb);
        
    end
    Ytilde(iii,:) = Ytilde_i;
    Mata_true(iii,:) = Mata_true_i;
    Matc_true(iii,:) = Matc_true_i;
    Matxi_true(iii) = Matxi_true_i;
end

% generate a new dataset using these household ids
boot_data = data;
for ii = 1:N
    
    % income
    boot_data(1+(ii-1)*7:ii*7,5) = Ytilde(ii,:);
    
    % consumption
    boot_data(1+(ii-1)*7:ii*7,9) = Matc_true(ii,:);
    
    % assets
    boot_data(1+(ii-1)*7:ii*7,12) = Mata_true(ii,:);
    
end

end