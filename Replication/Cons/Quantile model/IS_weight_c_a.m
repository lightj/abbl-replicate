function fval=IS_weight_c_a(i_Y,i_C,particle_draw,i_MatAGE_tot,MatzetaAGE1,ttt,...
                             Ntau,Vectau,Resqinit_eps,Resqinit,Resqinit_e0,Resqinit_c,...
                             b1_eps,bL_eps,b1,bL,b1_e0,bL_e0,b1_c,bL_c,Nis,T1,Mat,Mat2,...
                             Resqinit_a,b1_a,bL_a,Resqinit_a1,b1_a1,bL_a1,Mat3,A,asset_rule)

MatAGE_tot=i_MatAGE_tot(((ttt-1)*Nis+1:ttt*Nis),:);   
Matdraw=particle_draw;
Vect=i_Y(:,ttt)-Matdraw(:,ttt);

% Likelihood of the Data
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
dens4=zeros(Nis,1);

if ttt==T1
    
    
for jtau=1:Ntau-1
    dens2=dens2+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit_e0(:,jtau+1)-Resqinit_e0(:,jtau))'*MatzetaAGE1')'.*...
        (Matdraw(:,T1)>(Resqinit_e0(:,jtau)'*MatzetaAGE1')').*(Matdraw(:,T1)<=(Resqinit_e0(:,jtau+1)'*MatzetaAGE1')');
end

dens2=dens2+Vectau(1)*b1_e0*exp(b1_e0*(Matdraw(:,T1)-(Resqinit_e0(:,1)'*MatzetaAGE1')')).*...
    (Matdraw(:,T1)<=(Resqinit_e0(:,1)'*MatzetaAGE1')')+...
    (1-Vectau(Ntau))*bL_e0*exp(-bL_e0*(Matdraw(:,T1)-(Resqinit_e0(:,Ntau)'*MatzetaAGE1')')).*...
    (Matdraw(:,T1)>(Resqinit_e0(:,Ntau)'*MatzetaAGE1')');

if asset_rule==1
    for jtau=1:Ntau-1
        dens4=dens4+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit_a1(:,jtau+1)-Resqinit_a1(:,jtau))'*Mat3')'.*...
            (A(:,T1)>(Resqinit_a1(:,jtau)'*Mat3')').*(A(:,T1)<=(Resqinit_a1(:,jtau+1)'*Mat3')');
    end
    
    dens4=dens4+Vectau(1)*b1_a1*exp(b1_a1*(A(:,T1)-(Resqinit_a1(:,1)'*Mat3')')).*...
        (A(:,T1)<=(Resqinit_a1(:,1)'*Mat3')')+...
        (1-Vectau(Ntau))*bL_a1*exp(-bL_a1*(A(:,T1)-(Resqinit_a1(:,Ntau)'*Mat3')')).*...
        (A(:,T1)>(Resqinit_a1(:,Ntau)'*Mat3')');
end


% SUBSEQUENT PRIOR FOR PERIOD 2+
else
    

    for jtau=1:Ntau-1
        dens2=dens2+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit(:,jtau+1)-Resqinit(:,jtau))'*Mat')'.*...
            (Matdraw(:,ttt)>(Resqinit(:,jtau)'*Mat')').*(Matdraw(:,ttt)<=(Resqinit(:,jtau+1)'*Mat')');
    end
    
    dens2=dens2+Vectau(1)*b1*exp(b1*(Matdraw(:,ttt)-(Resqinit(:,1)'*Mat')')).*...
        (Matdraw(:,ttt)<=(Resqinit(:,1)'*Mat')')+...
        (1-Vectau(Ntau))*bL*exp(-bL*(Matdraw(:,ttt)-(Resqinit(:,Ntau)'*Mat')')).*...
        (Matdraw(:,ttt)>(Resqinit(:,Ntau)'*Mat')');
 
    if asset_rule==1
     for jtau=1:Ntau-1
        dens4=dens4+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit_a(:,jtau+1)-Resqinit_a(:,jtau))'*Mat3')'.*...
            (A(:,ttt)>(Resqinit_a(:,jtau)'*Mat3')').*(A(:,ttt)<=(Resqinit_a(:,jtau+1)'*Mat3')');
    end
    
    dens4=dens4+Vectau(1)*b1_a*exp(b1_a*(A(:,ttt)-(Resqinit_a(:,1)'*Mat3')')).*...
        (A(:,ttt)<=(Resqinit_a(:,1)'*Mat3')')+...
        (1-Vectau(Ntau))*bL_a*exp(-bL_a*(A(:,ttt)-(Resqinit_a(:,Ntau)'*Mat3')')).*...
        (A(:,ttt)>(Resqinit_a(:,Ntau)'*Mat3')');   
    end
    
end


% Likelihood of the Consumption data

Vect=i_C(:,ttt);


dens3=zeros(Nis,1);
for jtau=1:Ntau-1 
    dens3=dens3+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit_c(:,jtau+1)-Resqinit_c(:,jtau))'*Mat2')'.*...
        (Vect>(Resqinit_c(:,jtau)'*Mat2')').*(Vect<=(Resqinit_c(:,jtau+1)'*Mat2')');
end

temp =Vectau(1)*b1_c*exp(b1_c*(Vect-(Resqinit_c(:,1)'*Mat2')'));
temp(~isfinite(temp))=1000000000;
dens3=dens3+temp.*(Vect<=(Resqinit_c(:,1)'*Mat2')');

temp =(1-Vectau(Ntau))*bL_c*exp(-bL_c*(Vect-(Resqinit_c(:,Ntau)'*Mat2')'));
temp(~isfinite(temp))=1000000000;   
dens3=dens3+temp.*(Vect>(Resqinit_c(:,Ntau)'*Mat2')');  


if asset_rule==0
    dens4=ones(Nis,1);
end

fval=dens.*dens2.*dens3.*dens4;
end


