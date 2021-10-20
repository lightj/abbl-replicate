function Obj = wqregk_e0_age(c)
global tau Matdraw1 MatAGE1

Obj=mean((Matdraw1-MatAGE1*c).*...
    (tau-(Matdraw1-MatAGE1*c<0)));
    
end
