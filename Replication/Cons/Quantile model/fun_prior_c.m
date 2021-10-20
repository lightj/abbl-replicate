function fval=fun_prior_2factor(a_i_new,i_MatAGE1,Ntau,Vectau,Resqinit_zeta,b1_zeta,bL_zeta)

i_MatAGE1=i_MatAGE1(1,:);

denszeta=0;
for jtau=1:Ntau-1
    denszeta=denszeta+(Vectau(jtau+1)-Vectau(jtau))./((Resqinit_zeta(:,jtau+1)-Resqinit_zeta(:,jtau))'*i_MatAGE1')'.*...
        (a_i_new>(Resqinit_zeta(:,jtau)'*i_MatAGE1')').*(a_i_new<=(Resqinit_zeta(:,jtau+1)'*i_MatAGE1')');
end

denszeta=denszeta+Vectau(1)*b1_zeta*exp(b1_zeta*(a_i_new-(Resqinit_zeta(:,1)'*i_MatAGE1')')).*...
    (a_i_new<=(Resqinit_zeta(:,1)'*i_MatAGE1')')+...
    (1-Vectau(Ntau))*bL_zeta*exp(-bL_zeta*(a_i_new-(Resqinit_zeta(:,Ntau)'*i_MatAGE1')')).*...
    (a_i_new>(Resqinit_zeta(:,Ntau)'*i_MatAGE1')');

fval=denszeta;
end


