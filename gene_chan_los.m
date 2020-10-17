
xt = xt_s(:, t);
d2Dt = norm(xt);
d3Dt = sqrt(d2Dt^2+(-xt(3, 1)+hBS)^2);
Prlos = 1;
phit = atan(xt(2, 1)/xt(1, 1));
thetat = atan((-xt(3, 1)+hBS) / d2Dt) + pi * 0.5;
Dt = 4 * (hBS-1) * (xt(3, 1)-1) * fc * 1e9 / c;
if d2Dt <= Dt
    Plt = 32.4 + 21 * log10(d3Dt) + 20 * log10(fc);
else
    Plt = 32.4 + 40 * log10(d3Dt) + 20 * log10(fc) - 9.5 * log10(Dt^2 + (hBS - xt(3, 1))^2);
end
LOSindi = 1;
Pltp = 32.4 + 31.9 * log10(d3Dt) + 20 * log10(fc);
phiqt = phiqt_s(t, 1);
thetaqt = thetaqt_s(t, 1);
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
