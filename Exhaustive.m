%%Exhaustive
flag = 1;
snr_store = 0;
for p = 1 : Halfnumberb
    for pa = 1 : Halfnumberb
        for pppp = 1 : Halfnumberm
            for ppppa = 1 : Halfnumberm
                aMS = SMS(:, (pppp-1) * Halfnumberm + ppppa);
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
    end
end
if flag < 0.5
    ahievrate(3, nnids) = ahievrate(3, nnids) + log2(1 + snr_store);
    succcount(3,1) = succcount(3,1) + 1;
else
    misdetection(3, nnids) = misdetection(3, nnids) + 1;
end