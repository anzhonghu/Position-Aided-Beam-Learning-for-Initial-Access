%%%PropSVM
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
if t<=Tsetnumber
    for p = 1 : Halfnumberb
        for pa = 1 : Halfnumberb
            for pppp = 1 : Halfnumberm
                for ppppa = 1 : Halfnumberm
                    aMS = SMS(:, (pppp-1) * Halfnumberm + ppppa);
                    aBS = SBS(:, (p-1) * Halfnumberb + pa);
                    beamstore_all = (p-1) * Halfnumberb + pa;
                    yt_s = abs(aBS' * Hmultipath * aMS)^2 * 10^(MSpower*0.1) /  (10^(0.1*Noisepower)*abs(randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5)^2);
                    if yt_s > 10^(0.1*Thetathresdb)
                        flag = 0;
                        if snr_store < yt_s
                            snr_store = yt_s;
                            MSbeamupdate_record = aMS;
                        else
                        end
                    else
                    end
                end
            end
        end
    end
    xtsto(:, t) = xt;
    if flag > 0.5
        yqqa(:,t) = -1 * ones(numberb,1);
    else
        aMS = MSbeamupdate_record;
        for p = 1 : Halfnumberb
            for pa = 1 : Halfnumberb
                aBS = SBS(:, (p-1) * Halfnumberb + pa);
                beamstore_all = (p-1) * Halfnumberb + pa;
                yt_s = abs(aBS' * Hmultipath * aMS)^2 * 10^(MSpower*0.1) /  (10^(0.1*Noisepower)*abs(randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5)^2);
                if yt_s > 10^(0.1*Thetathresdb)
                    yqqa(beamstore_all,t) = 1;
                else
                    yqqa(beamstore_all,t) = -1;
                end
            end
        end
    end
    %     if flag < 0.5%%success
    %         access_delay(2, nnids) = access_delay(2, nnids) + numberm * numberb * Tper + numberb * Tper;
    %         succcount(5,1) = succcount(5,1) + 1;
    %     else%%fail
    %         misdetection(5, nnids) = misdetection(5, nnids) + 1;
    %     end
else
    beam_fxxts = zeros(numberb, 2);
    for paaa = 1 : Halfnumberb
        for paaaa = 1 : Halfnumberb
            fxxt = 0;
            beamstore_all = (paaa-1) * Halfnumberb + paaaa;
            for txx = 1 : Tsetnumber
                fxxt = fxxt + yqqa(beamstore_all,txx) * alphaqqa(beamstore_all,txx) * exp(-norm(xtsto(:, txx)-xt)^2/sigmasvm);
            end
            beam_fxxts(beamstore_all, :) = [beamstore_all fxxt];
        end
    end
    [~, orderind] = sort(beam_fxxts(:,2),'descend');
    for s = 1 : SS
        fxxt = beam_fxxts(orderind(s,1), 2);
        beamstore_all = beam_fxxts(orderind(s,1), 1);
        if fxxt > 0
            %%%%%%%%%%%%%%%%%%%%%%%%
            aBS = SBS(:, beamstore_all);
            for p = 1 : Halfnumberm
                for pa = 1 : Halfnumberm
                    aMS = SMS(:, (p-1) * Halfnumberm + pa);
                    yt_s = abs(aBS' * Hmultipath * aMS)^2 * 10^(MSpower*0.1) /  (10^(0.1*Noisepower)*abs(randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5)^2);
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
        else
        end
    end
    if flag < 0.5%%success
        access_delay(2, nnids) = access_delay(2, nnids) + numberm * SS * Tper + SS * Tper;
        succcount(5,1) = succcount(5,1) + 1;
        ahievrate(5, nnids) = ahievrate(5, nnids) + log2(1 + snr_store);
    else%%fail
        misdetection(5, nnids) = misdetection(5, nnids) + 1;
    end
end
if t==Tsetnumber
    Kcor = zeros(Tsetnumber, Tsetnumber);
    for t11 = 1 : Tsetnumber
        for t22 = 1 : Tsetnumber
            Kcor(t11, t22) = exp(-norm(xtsto(:, t11)-xtsto(:, t22))^2/sigmasvm);
        end
    end
    for p = 1 : Halfnumberb
        for pa = 1 : Halfnumberb
            bsbeamcountforopti = (p-1) * Halfnumberb + pa;
            Dqqa = (yqqa(bsbeamcountforopti, 1:Tsetnumber)' * yqqa(bsbeamcountforopti, 1:Tsetnumber)) .* Kcor;
            %%%%%%%%%%%%%%%%%%%%%%%
            Dqqaqq = [Dqqa zeros(Tsetnumber,1); zeros(Tsetnumber,1)' 1/C];
            fsolv = [-ones(Tsetnumber,1); 0]';
            Asolv = [eye(Tsetnumber) -ones(Tsetnumber,1);zeros(1,Tsetnumber+1)];
            bsolv = zeros(Tsetnumber+1,1);
            Aeqsolv = [yqqa(bsbeamcountforopti, 1:Tsetnumber) 0];
            beqsolv = 0;
            lbsolv = zeros(Tsetnumber+1,1);
            xsolver = quadprog(Dqqaqq,fsolv,Asolv,bsolv,Aeqsolv,beqsolv,lbsolv,[]);
            alphaqqa(bsbeamcountforopti,1:Tsetnumber) = xsolver(1:Tsetnumber,1).';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end


