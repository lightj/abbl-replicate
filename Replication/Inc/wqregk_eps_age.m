function Obj = wqregk_eps_age(c)
global tau Matdraw_tot MatAGE_tot Y_tot Matdraw_eps

Obj=mean((Matdraw_eps-MatAGE_tot*c).*...
    (tau-(Matdraw_eps-MatAGE_tot*c<0)));
    
end
