function Obj = wqregk_cons_age(c)
global tau C_tot MatClag

Obj=mean((C_tot-MatClag*c).*...
    (tau-(C_tot-MatClag*c<0)));
    
end
