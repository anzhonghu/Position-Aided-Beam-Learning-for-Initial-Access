%%%CIin3
pccc = round((phiqt+pi*0.5) * Halfnumberm / pi+1);
if pccc<1
    pccc = 1;
end
if pccc > Halfnumberm
    pccc = Halfnumberm;
end
paccc = round(thetaqt * Halfnumberm / pi + 1);
if paccc<1
    paccc = 1;
end
if paccc > Halfnumberm
    paccc = Halfnumberm;
end
posi_in_potenbeam = (pccc-1) * Halfnumberm + paccc;

aMS = SMS(:, posi_in_potenbeam);
flag = 1;
snr_store = 0;
for p = 1 : Halfnumberb
    for pa = 1 : Halfnumberb
        aBS = SBS(:, (p-1) * Halfnumberb + pa);
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
if flag < 0.5
    ahievrate(1, nnids) = ahievrate(1, nnids) + log2(1 + snr_store);
    succcount(1,1) = succcount(1,1) + 1;
else
    misdetection(1, nnids) = misdetection(1, nnids) + 1;
end