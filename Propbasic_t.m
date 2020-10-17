%%%Propbasic
pccc = round((phit+pi/3) * Halfnumberb / pi / 2 * 3+1);
if pccc<1
    pccc = 1;
end
if pccc > Halfnumberb
    pccc = Halfnumberb;
end
paccc = round(thetat * Halfnumberb / pi + 1);
if paccc<1
    paccc = 1;
end
if paccc > Halfnumberb
    paccc = Halfnumberb;
end
posi_in_potenbeam = (pccc-1) * Halfnumberb + paccc;

flag = 1;
snr_store = 0;

if t <=Tsetnumber
    for p = 1 : Halfnumberb
        for pa = 1 : Halfnumberb
            for pppp = 1 : Halfnumberm
                for ppppa = 1 : Halfnumberm
                    aMSx = SMS(:, (pppp-1) * Halfnumberm + ppppa);
                    aBSx = SBS(:, (p-1) * Halfnumberb + pa);
                    if LOS_indi > 0.5
                        sigmatdb = sigmatdb_s(t, 1);
                        betat = 10^(0.1*sigmatdb) * (randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5);
                        Hmultipath = betat * aBS * aMS';
                    else
                    end
                    yt_s = abs(aBSx' * Hmultipath * aMSx)^2 * 10^(MSpower*0.1) / (10^(0.1*Noisepower)*abs(randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5)^2);
                    phipacz = -pi/3 + (p-1) / Halfnumberb * pi * 2 / 3;
                    thetapacz = (pa-1) / Halfnumberb * pi;
                    disz = abs(sin(phit)*sin(thetat)-sin(phipacz)*sin(thetapacz))+abs(cos(thetat)-cos(thetapacz));
                    if yt_s > 10^(0.1*Thetathresdb) 
                        flag = 0;
                        if snr_store < yt_s
                            snr_store = yt_s;
                            beamstore_fail = (p-1) * Halfnumberb + pa;
                        else
                        end
                    else
                    end
                end
            end
        end
    end
    if flag < 0.5%%success
        if potenbeamcount(posi_in_potenbeam,1) < S
            potenbeamcount(posi_in_potenbeam,1) = potenbeamcount(posi_in_potenbeam,1) + 1;
            SBS_Basic(:, (posi_in_potenbeam-1)*S+potenbeamcount(posi_in_potenbeam,1)) = SBS(:, beamstore_fail);
            SBS_Basic_indi(1, (posi_in_potenbeam-1)*S+potenbeamcount(posi_in_potenbeam,1)) = beamstore_fail;
        else
        end
    else
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    for s = 1 : S
        if SBS_Basic_indi(1, (posi_in_potenbeam-1)*S+s) > 0
            aBSx = SBS_Basic(:, (posi_in_potenbeam-1)*S+s);
        else
        end
        for p = 1 : Halfnumberm
            for pa = 1 : Halfnumberm
                aMSx = SMS(:, (p-1) * Halfnumberm + pa);
                if LOS_indi > 0.5
                    sigmatdb = sigmatdb_s(t, 1);
                    betat = 10^(0.1*sigmatdb) * (randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5);
                    Hmultipath = betat * aBS * aMS';
                else
                end
                yt_s = abs(aBSx' * Hmultipath * aMSx)^2 * 10^(MSpower*0.1) /  (10^(0.1*Noisepower)*abs(randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5)^2);
                if yt_s > 10^(0.1*Thetathresdb)
                    flag = 0;
                    if snr_store < yt_s
                        snr_store = yt_s;
                    else
                    end
                else
                end
            end
        end
    end
    if flag < 0.5%%success
        ahievrate(4, nnids) = ahievrate(4, nnids) + log2(1 + snr_store);
    else%%fail
        misdetection(4, nnids) = misdetection(4, nnids) + 1;
    end
end
