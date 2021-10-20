function Obj = wqregk_a_age(c)
global tau A1 MatA

Obj=mean((A1-MatA*c).*...
    (tau-(A1-MatA*c<0)));
    
end
