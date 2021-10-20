function fval=IS_weight(i_Y,particle_draw,i_MatAGE_t,MatAGE1,AGE,ttt,...
                             Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,...
                             b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,K1,K2,meanAGE,...
                             stdAGE,meanY,stdY,Nis,T1,Mat)

MatAGE_tot=i_MatAGE_t(((ttt-1)*Nis+1:ttt*Nis),:);   
Matdraw=particle_draw;
Vect=i_Y(:,ttt)-Matdraw(:,ttt);

%Likelihood of the Data
dens=zeros(Nis,1);
for jtau=1:Ntau-1 
    dens=dens+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit_eps(:,jtau+1)-Resqinit_eps(:,jtau))'*MatAGE_tot')'.*...
        (Vect>(Resqinit_eps(:,jtau)'*MatAGE_tot')').*(Vect<=(Resqinit_eps(:,jtau+1)'*MatAGE_tot')');
end

dens=dens+Vectau(1)*b1_eps*exp(b1_eps*(Vect-(Resqinit_eps(:,1)'*MatAGE_tot')')).*...
    (Vect<=(Resqinit_eps(:,1)'*MatAGE_tot')')+...
    (1-Vectau(Ntau))*bL_eps*exp(-bL_eps*(Vect-(Resqinit_eps(:,Ntau)'*MatAGE_tot')')).*...
    (Vect>(Resqinit_eps(:,Ntau)'*MatAGE_tot')');


% PRIOR FOR PERIOD 1
dens2=zeros(Nis,1);
if ttt==T1
    
for jtau=1:Ntau-1
    dens2=dens2+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit_e0(:,jtau+1)-Resqinit_e0(:,jtau))'*MatAGE1')'.*...
        (Matdraw(:,T1)>(Resqinit_e0(:,jtau)'*MatAGE1')').*(Matdraw(:,T1)<=(Resqinit_e0(:,jtau+1)'*MatAGE1')');
end

dens2=dens2+Vectau(1)*b1_e0*exp(b1_e0*(Matdraw(:,T1)-(Resqinit_e0(:,1)'*MatAGE1')')).*...
    (Matdraw(:,T1)<=(Resqinit_e0(:,1)'*MatAGE1')')+...
    (1-Vectau(Ntau))*bL_e0*exp(-bL_e0*(Matdraw(:,T1)-(Resqinit_e0(:,Ntau)'*MatAGE1')')).*...
    (Matdraw(:,T1)>(Resqinit_e0(:,Ntau)'*MatAGE1')');


% SUBSEQUENT PRIOR FOR PERIOD 2+
else
    
tt=ttt-1;

    for jtau=1:Ntau-1
        dens2=dens2+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit(:,jtau+1)-Resqinit(:,jtau))'*Mat')'.*...
            (Matdraw(:,tt+1)>(Resqinit(:,jtau)'*Mat')').*(Matdraw(:,tt+1)<=(Resqinit(:,jtau+1)'*Mat')');
    end
    
    dens2=dens2+Vectau(1)*b1*exp(b1*(Matdraw(:,tt+1)-(Resqinit(:,1)'*Mat')')).*...
        (Matdraw(:,tt+1)<=(Resqinit(:,1)'*Mat')')+...
        (1-Vectau(Ntau))*bL*exp(-bL*(Matdraw(:,tt+1)-(Resqinit(:,Ntau)'*Mat')')).*...
        (Matdraw(:,tt+1)>(Resqinit(:,Ntau)'*Mat')');
    
end


fval=dens.*dens2;
end


