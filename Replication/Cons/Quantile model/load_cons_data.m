function load_cons_data(data)

global T N K1 K2 K3 K3t K4 K5 Ntau Vectau tau ...
       D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
       Matdraw1 AGE1 Y1 MatAGE1 Matdraw_t Matdraw_lag...
       AGE_t Matdraw_t_lag Y_tot Matdraw_tot MatAGE_tot Y_t Y_lag...
       AGE_tot meanT stdT C A A_tot  C_tot meanA stdA meanC stdC A_t A_lag C_lag ...
       Vect_c_xi  Vect_c_age Vect_c_eta Vect_c_eps Vect_c_a  A1 Vect_a_xi Vect_a_age Vect_a_eta ...
       Vect_a_eps Vect_a_c Vect_a_a Vect_xi_a Vect_xi_age Vect_xi_eta Vect_a1_eta Vect_a1_age...
        Knots_c_xi Knots_c_age Knots_c_eta Knots_c_eps Knots_c_a Knots_a_xi Knots_a_age Knots_a_eta ...
        Knots_a_eps Knots_a_c Knots_a_a Knots_xi_a Knots_xi_age Knots_xi_eta Knots_a1_eta Knots_a1_age...
        M12 MatXi M12b EDUC YB meanEDUC meanYB stdEDUC stdYB Ybar M12a meanYbar stdYbar
    
C=data(:,9);
MatC=zeros(N,T);
for tt=1:T
    MatC(:,tt)=C(tt:T:N*T);
end
C=MatC;

A=data(:,12);
MatA=zeros(N,T);
for tt=1:T
    MatA(:,tt)=A(tt:T:N*T);
end
A=MatA;

EDUC=data(:,6);
MatEDUC=zeros(N,T);
for tt=1:T
    MatEDUC(:,tt)=EDUC(tt:T:N*T);
end

YB=data(:,7);
MatYB=zeros(N,T);
for tt=1:T
    MatYB(:,tt)=YB(tt:T:N*T);
end
YB=zeros(N,1);

EDUC=zeros(N,1);
for iii=1:N
    EDUC(iii)=MatEDUC(iii,ENT(iii));
    YB(iii)=MatYB(iii,ENT(iii));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    _TOT (for consumption)                                    %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample=find(D(:));

A_tot=A(:);
A_tot=A_tot(sample);
meanA=mean(A_tot);
stdA=std(A_tot);

C_tot=C(:);
C_tot=C_tot(sample);
meanC=mean(C_tot);
stdC=std(C_tot);

AGE_tot=AGE(:);
AGE_tot=AGE_tot(sample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      _T  (for assets)     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample=find(D_t(:));

A_t=A(:,2:T);
A_t=A_t(:);
A_t=A_t(sample);


A_lag=A(:,1:T-1);
A_lag=A_lag(:);
A_lag=A_lag(sample);

C_lag=C(:,1:T-1);
C_lag=C_lag(:);
C_lag=C_lag(sample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      _1  (for assets)     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = zeros(N,1);
AGE1 = zeros(N,1);
Y1 = zeros(N,1);
Ybar =zeros(N,1);
for iii=1:N
    Y1(iii)=Y(iii,ENT(iii));
    A1(iii)=A(iii,ENT(iii));
    AGE1(iii)=AGE(iii,ENT(iii));
    Ybar(iii)=mean(Y(iii,ENT(iii):EXT(iii)));
end

%bols = pinv([ones(N,1) AGE1 AGE1.^2])*Ybar;
%Ybar =Ybar - ( [ones(N,1) AGE1 AGE1.^2]*bols);

meanYbar =mean(Ybar);
stdYbar=std(Ybar);
meanEDUC=mean(EDUC);
stdEDUC=std(EDUC);
meanYB=mean(YB);
stdYB=std(YB);

MatXi=[];
for kk12=0:M12
    for kk13=0:M12a
        for kk14=0:M12b
        MatXi=[MatXi hermite(kk12,(EDUC-meanEDUC)/stdEDUC)...
            .*hermite(kk13,(YB-meanYB)/stdYB)...
            .*hermite(kk14,(Ybar-meanYbar)/stdYbar)];
        end
    end
end

end