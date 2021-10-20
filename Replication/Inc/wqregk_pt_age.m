function Obj = wqregk_pt_age(c)
global tau Matdraw_t Mat_ETA_AGE 

Obj=mean((Matdraw_t-Mat_ETA_AGE*c).*...
    (tau-(Matdraw_t-Mat_ETA_AGE*c<0)));
    
end
