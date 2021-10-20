function Obj = wqregk_xi_age(c)
global tau Matdraw MatXi T

Obj=mean((Matdraw(:,T+1)-MatXi*c).*...
    (tau-(Matdraw(:,T+1)-MatXi*c<0)));
    
end
