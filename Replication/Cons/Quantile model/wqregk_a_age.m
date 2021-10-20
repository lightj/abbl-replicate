function Obj = wqregk_a_age(c)
global tau A_t MatAlag

Obj=mean((A_t-MatAlag*c).*...
    (tau-(A_t-MatAlag*c<0)));
    
end
