function load_data(data)

global T N K1 K2 K3 K3t K4 K5 K6 Ntau Vectau tau ...
       D D_t Y AGE meanAGE stdAGE meanY stdY ENT EXT...
       Matdraw1 AGE1 Y1 MatAGE1 Matdraw_t Matdraw_lag...
       AGE_t Matdraw_t_lag Y_tot Matdraw_tot MatAGE_tot Y_t Y_lag...
       AGE_tot meanT stdT  EDUC YB meanEDUC meanYB stdEDUC stdYB


Y=data(:,5);

N=size(Y,1)/T;
MatY=zeros(N,T);
for tt=1:T
    MatY(:,tt)=Y(tt:T:N*T);
end
Y=MatY;

D=data(:,4);
MatD=zeros(N,T);
for tt=1:T
    MatD(:,tt)=D(tt:T:N*T);
end
D=MatD;

D_t=D(:,1).*D(:,2);
for j=3:T
    D_t=[D_t;D(:,j-1).* D(:,j)];
end

AGE=data(:,3);
MatAGE=zeros(N,T);
for tt=1:T
    MatAGE(:,tt)=AGE(tt:T:N*T);
end
AGE=MatAGE;

[~,ENT]=max(D>0,[],2);
EXT=sum(D,2);
EXT=ENT+EXT-1;

% meanT=mean(ENT);
% stdT=std(ENT);

meanT=0;
stdT=1;

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

meanEDUC=mean(EDUC);
stdEDUC=std(EDUC);
meanYB=mean(YB);
stdYB=std(YB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    _TOT (for epsilon)   %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample=find(D(:));

AGE_tot=AGE(:);
AGE_tot=AGE_tot(sample);
meanAGE=mean(AGE_tot);
stdAGE=std(AGE_tot);

Time=[];
for tt=1:T
 Time=[Time tt*ones(N,1)]; 
end
Time_tot=Time(:);
Time_tot=Time_tot(sample);

Y_tot=Y(:);
Y_tot=Y_tot(sample);
meanY=mean(Y_tot);
stdY=std(Y_tot);

MatAGE_tot=[];
for kk3=0:K3
    for kk3t=0:K3t
    MatAGE_tot=[MatAGE_tot hermite(kk3,(AGE_tot-meanAGE)/stdAGE)...
        .*hermite(kk3t,(Time_tot-meanT)/stdT)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      _1 (for eta1)      %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AGE1=zeros(N,1);
Y1=zeros(N,1);
for ii=1:N
AGE1(ii)=AGE(ii,ENT(ii));
Y1(ii)=Y(ii,ENT(ii));    
end

MatAGE1=[];
for kk4=0:K4
    for kk5=0:K5
        for kk6=0:K6
    MatAGE1=[MatAGE1 hermite(kk4,(AGE1-meanAGE)/stdAGE)...
        .*hermite(kk5,(EDUC-meanEDUC)/stdEDUC)...
        .*hermite(kk6,(YB-meanYB)/stdYB)];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      _T  (for eta)     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample=find(D_t(:));

AGE_t=AGE(:,2);
for j=3:T
    AGE_t=[AGE_t;AGE(:,j)];
end
AGE_t=AGE_t(sample);

Y_t=Y(:,2);
for j=3:T
    Y_t=[Y_t;Y(:,j)];
end
Y_t=Y_t(sample);

Y_lag=Y(:,1:T-1);
Y_lag=Y_lag(:);
Y_lag=Y_lag(sample);
end