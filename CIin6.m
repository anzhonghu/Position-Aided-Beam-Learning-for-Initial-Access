%%%CIin6
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
aBS = SBS_Basic(:, (posi_in_potenbeam-1)*S+1);

flag = 1;
snr_store = 0;
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
if flag < 0.5
    ahievrate(2, nnids) = ahievrate(2, nnids) + log2(1 + snr_store);
    succcount(2,1) = succcount(2,1) + 1;
else
    misdetection(2, nnids) = misdetection(2, nnids) + 1;
end