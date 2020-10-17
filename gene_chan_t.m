
xt = zeros(3, 1);
xt(1:2, 1) = Rmax * 2 * (rand(2, 1)-0.5*ones(2,1));
while norm(xt) > Rmax || norm(xt) < Rmin || xt(1,1)<0 || abs(atan(xt(2,1)/xt(1,1)))>pi/3
    xt(1:2, 1) = Rmax * 2 * (rand(2, 1)-0.5*ones(2,1));
end
d2Dt = norm(xt);
xt(3, 1) = 1.5 + (22.5-1.5) * rand(1, 1);
d3Dt = norm(xt);
if d2Dt > 18
    Prlos = 18 / d2Dt + (1-18 / d2Dt) * exp(-d2Dt/36);
else
    Prlos = 1;
end
phit = atan(xt(2, 1)/xt(1, 1));
thetat = atan((-xt(3, 1)+hBS) / d2Dt) + pi * 0.5;
Lostv = rand(1, 1);
Dt = 4 * (hBS-1) * (xt(3, 1)-1) * fc * 1e9 / c;
if d2Dt <= Dt
    Plt = 32.4 + 21 * log10(d3Dt) + 20 * log10(fc);
else
    Plt = 32.4 + 40 * log10(d3Dt) + 20 * log10(fc) - 9.5 * log10(Dt^2 + (hBS - xt(3, 1))^2);
end
SFt = randn(1, 1) * sigmasfdb;
%     Adbv = -min(12 * ((thetat*180/pi-90)/theta3dbdegree)^2, SLAVdb);
%     Adbh = -min(12 * (phit*180/pi/phi3dbdegree)^2, Amaxdb);
%     Adb = GEmax-min(-(Adbv+Adbh), Amaxdb);
sigmatdb = -Plt - SFt;% + Adb;
if Lostv <= Prlos
    LOSindi = 1;
else
    LOSindi = 0;
end
Pltp = 32.4 + 31.9 * log10(d3Dt) + 20 * log10(fc);
phiqt = rand(1,1) * pi;
thetaqt = rand(1,1) * pi * 0.5;
Hmultipath = zeros(Nbsh*Nbsv, Nmsh*Nmsv);
sigmanclnraysumfpmd = zeros(Ncl, 1);
aMSs = zeros(Nmsh*Nmsv, Ncl*Nray+1);
for ncl = 1 : Ncl
    Ynclt = randn(1,1) * ASAdegree / 7 * pi / 180;
    SFtnclp =  randn(1, 1) * sigmasfadb;
    if xt(1,1)<=Rmid && xt(2,1)>=0
        phitclt = clangaz(1, ncl);
        thetatnclt = clangel(1, ncl);
        phiqtclt = clangaz(1, ncl);
        thetaqtclt = clangel(1, ncl);
    end
    if xt(1,1)<=Rmid && xt(2,1)<0
        phitclt = clangaz(2, ncl);
        thetatnclt = clangel(2, ncl);
        phiqtclt = clangaz(2, ncl);
        thetaqtclt = clangel(2, ncl);
    end
    if xt(1,1)>Rmid && xt(2,1)>=0
        phitclt = clangaz(3, ncl);
        thetatnclt = clangel(3, ncl);
        phiqtclt = clangaz(3, ncl);
        thetaqtclt = clangel(3, ncl);
    end
    if xt(1,1)>Rmid && xt(2,1)<0
        phitclt = clangaz(4, ncl);
        thetatnclt = clangel(4, ncl);
        phiqtclt = clangaz(4, ncl);
        thetaqtclt = clangel(4, ncl);
    end
    for nray = 1 : Nray
        phiqtnclnray = phiqtclt + cASAdegree * alphatabledegree(nray, 1) * pi / 180;
        thetaqtnclnray = thetaqtclt + cASAdegree * alphatabledegree(nray, 1) * pi / 180;
        phitnclnray = phitclt + cASAdegree * alphatabledegree(nray, 1) * pi / 180;
        thetatnclnray = thetatnclt + cASAdegree * alphatabledegree(nray, 1) * pi / 180;
        %         Adbv = -min(12 * ((thetatnclnray*180/pi-90)/theta3dbdegree)^2, SLAVdb);
        %         Adbh = -min(12 * (phitnclnray*180/pi/phi3dbdegree)^2, Amaxdb);
        %         Adb = GEmax -min(-(Adbv+Adbh), Amaxdb);
        sigmanclnray = -Pltp - SFtnclp;% + Adb;
        betanclnray = (randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5) * 10^(0.1*sigmanclnray);
        sigmanclnraysumfpmd(ncl, 1) = 10^(0.2*sigmanclnray);
        for nbshcc = 1 : Nbsh
            for nbsvcc = 1 : Nbsv
                aBS((nbshcc-1)*Nbsh+nbsvcc, 1) = 1 / sqrt(Nbsh*Nbsv)...
                    *exp(1i * dvlambda * 2 * pi * ((nbsvcc-1)*sin(phitnclnray)*sin(thetatnclnray)+(nbshcc-1)*cos(thetatnclnray)));
            end
        end
        for nmshcc = 1 : Nmsh
            for nmsvcc = 1 : Nmsv
                aMS((nmshcc-1)*Nmsh+nmsvcc, 1) = 1 / sqrt(Nmsh*Nmsv)...
                    *exp(1i * dvlambda * 2 * pi * ((nmsvcc-1)*sin(phiqtnclnray)*sin(thetaqtnclnray)+(nmshcc-1)*cos(thetaqtnclnray)));
            end
        end
        aMSs(:, (ncl-1)*Nray+nray) = aMS;
        Hmultipath = Hmultipath + betanclnray * aBS * aMS';
    end
end
if LOSindi > 0.5
    for nbshcc = 1 : Nbsh
        for nbsvcc = 1 : Nbsv
            aBS((nbshcc-1)*Nbsh+nbsvcc, 1) = 1 / sqrt(Nbsh*Nbsv)...
                *exp(1i * dvlambda * 2 * pi * ((nbsvcc-1)*sin(phit)*sin(thetat)+(nbshcc-1)*cos(thetat)));
        end
    end
    for nmshcc = 1 : Nmsh
        for nmsvcc = 1 : Nmsv
            aMS((nmshcc-1)*Nmsh+nmsvcc, 1) = 1 / sqrt(Nmsh*Nmsv)...
                *exp(1i * dvlambda * 2 * pi * ((nmsvcc-1)*sin(phiqt)*sin(thetaqt)+(nmshcc-1)*cos(thetaqt)));
        end
    end
    aMSs(:,end) = aMS;
    betat = 10^(0.1*sigmatdb) * (randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5);
    Hmultipath = sqrt(1 / (Kr+1)) * betat * aBS * aMS' + sqrt(Kr / (Kr+1)) * Hmultipath;
else
end