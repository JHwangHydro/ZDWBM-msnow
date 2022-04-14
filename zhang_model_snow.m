function[Qest, Evap, Storage, Sf, Gf, SSf, Snowpack] = zhang_model_snow(P,PET,SF,nmonths,beta,alpha1,alpha2,Smax,d,So,Go,SSo, Tmax, Tmin)

Qest = zeros(nmonths,1);
Evap = zeros(nmonths,1);
Storage = zeros(nmonths,1);
Snowpack = zeros(nmonths,1);

SSini = SSo;
Sini = So;
Gini = Go;

for i = 1:nmonths

    k = rem(i,12);
    if k == 0
        k = 12;
    end
    
    RFt = P(i);
    SFt = SF(i);
    Tmaxt = Tmax(i);
    Tmint = Tmin(i);
    
    if RFt + SFt < 0.01
        RFt = 0.01;
    end
    
%   St-1 Freezing
    Tavg = (Tmaxt + Tmint) * 0.5 * 0.1;
    if Tavg < -0.5
        SSini = SSini + Sini;
        Sini = 0;
    end
    
%   Radiative Heat Flux Seperated for PET and PM    
    a = (RFt + Sini) / (RFt + SFt + Sini + SSini);
    b = (SFt + SSini) / (RFt + SFt + Sini + SSini);
    
    if PET(i) == 0
       PET(i) = 0.1;
    end
    
    PETt = PET(i) * a;
    PMt = PET(i) * (b/0.15);
    
    Mt = (SFt + SSini)*(1 + (PMt/(SFt+SSini)) - (1 + (PMt/(SFt + SSini)).^(1/(1-beta(k)))).^(1-beta(k)));
    
    if isnan(Mt)
        Mt = 0;
    end
    
    SSt = SSini + SFt - Mt;
    if SSt < 0
        SSt = 0;
    end
    
    melt(i) = Mt/(SFt + SSini);
    pmelt(i) = PMt/(SFt + SSini);
    if isnan(melt(i))
        melt(i) = 0;
    end
    if isnan(pmelt(i))
        pmelt(i) = 0;
    end
    
    Snowpack(i) = SSt;
    Pt = RFt + Mt;
    
    Xot = Smax(k) - Sini + PETt;

    if Xot < 0
        Xot = 0;
    end
    
    if Pt < 0.01
        Pt = 0.01;
    end
    
    Xt = Pt*(1+ (Xot/Pt) - (1 + (Xot/Pt).^(1/(1-alpha1(k)))).^(1-alpha1(k)));

    Wt = Xt + Sini;
    
    if Wt == 0
        Wt = 0.001;
    end

    Yt = Wt*(1 + ((Smax(k)+PETt)/Wt) - (1 + ((Smax(k)+PETt)/Wt).^(1/(1-alpha2(k)))).^(1-alpha2(k)));

    Et = Wt*(1 + (PETt/Wt) - (1 + (PETt/Wt).^(1/(1-alpha2(k)))).^(1-alpha2(k)));

    St = Yt - Et;
    
    if St < 0
        St = 0;
    end    

    Rt = Wt - Yt;

    Gt = (1-d(k))*Gini + Rt;

    runoff_direct = Pt - Xt;

    runoff_base = d(k)*Gini;

    Qest(i) = runoff_direct + runoff_base;
    Evap(i) = Et;
    Storage(i) = St;

    
    Sini = St;
    Gini = Gt;
    SSini = SSt;
    
    Pre(i,1) = Pt;
    Pevap(i,1) = PETt;
    Ret(i,1) = Xt;
    Ret_lim(i,1) = Xot;
    Wavail(i,1) = Wt;
    Evapop(i,1) = Yt;
    Evapor(i,1) = Et;
    Stor(i,1) = St;
    Rechar(i,1) = Rt;
    Ground(i,1) = Gt;
    ro_dir(i,1) = runoff_direct;
    ro_base(i,1) = runoff_base; 
    
    Sf = St;
    Gf = Gt;
    SSf = SSt;
end

return