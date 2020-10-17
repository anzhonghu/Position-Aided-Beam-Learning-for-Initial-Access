close all;
clear;
runId = 1;
rng(runId);
LOS_indi = 1;
T = 10;
S = 4;
SS = 4;
Amaxdb = 30;
theta3dbdegree = 65;
SLAVdb = 30;
phi3dbdegree = 65;
GEmax = 8;
dvlambda = 0.5;
hBS = 10;
Kr = 10^0.9;
Tper = 5;
fc = 28;
c = 3*1e8;
sigmasfdb = 4;
sigmasfadb = 8.2;
ASAdegree = 14.8;
cASAdegree = 22;
Nbsh = 8;
Nbsv = 8;
Nmsh = 2;
Nmsv = 2;
aBS = zeros(Nbsh*Nbsv, 1);
aMS = zeros(Nmsh*Nmsv, 1);
Bigthetadb = -4;
C = 1;
sigmasvm = 1e2;
Ttotal = 1e4;
Rmin = 10;
Rmax = 100;
Rmid = (Rmax+Rmin) * 0.5;
Thetathresdb = -4;
Noisepower = -174;
alphatabledegree = [0.0447;-0.0447;0.1413;-0.1413;0.2492;-0.2492;0.3715;-0.3715;...
    0.5129;-0.5129;0.6797;-0.6797;0.8844;-0.8844;1.1481;-1.1481; 1.5195;-1.5195;...
    2.1551;-2.1551];
%divide into four areas
%clusterangleazimuth, 4 rows
clangaz = [-60 -59 -58 -57 -56; -30 -29 -28 -27 -26; 30 29 28 27 26; 60 59 58 57 56];
%clusterangleelevation, 4 rows
clangel = [-60 -59 -58 -57 -56; -30 -29 -28 -27 -26; 30 29 28 27 26; 60 59 58 57 56];
BSpower = 23;
MSpower = 23;
Qm = 4;
Halfnumberm = 2 ^(0.5*Qm);
numberm = 2 ^Qm;
SMS = zeros(Nmsh*Nmsv, 2^Qm);
for p = 1 : Halfnumberm
    for pa = 1 : Halfnumberm
        phipac = -pi*0.5 + (p-1) / Halfnumberm * pi;
        thetapac = (pa-1) / Halfnumberm * pi;
        for nmshcc = 1 : Nmsh
            for nmsvcc = 1 : Nmsv
                aMS((nmshcc-1)*Nmsh+nmsvcc, 1) = 1 / sqrt(Nmsh*Nmsv)...
                    *exp(1i * dvlambda * 2 * pi * ((nmsvcc-1)*sin(phipac)*sin(thetapac)+(nmshcc-1)*cos(thetapac)));
            end
        end
        SMS(:, (p-1)*Halfnumberm + pa) = aMS;
    end
end
Ss = [2;4;6;8;10];
Qb = 6;
Halfnumberb = 2 ^(0.5*Qb);
numberb = 2 ^Qb;
SBS = zeros(Nbsh*Nbsv, 2^Qb);
for p = 1 : Halfnumberb
    for pa = 1 : Halfnumberb
        phipac = -pi/3 + (p-1) / Halfnumberb * pi * 2 / 3;
        thetapac = (pa-1) / Halfnumberb * pi;
        for nbshcc = 1 : Nbsh
            for nbsvcc = 1 : Nbsv
                aBS((nbshcc-1)*Nbsh+nbsvcc, 1) = 1 / sqrt(Nbsh*Nbsv)...
                    *exp(1i * dvlambda * 2 * pi * ((nbsvcc-1)*sin(phipac)*sin(thetapac)+(nbshcc-1)*cos(thetapac)));
            end
        end
        SBS(:, (p-1) * Halfnumberb + pa) = aBS;
    end
end
misdetection = zeros(5, length(Ss));%5 methods
Tsetnumber = 100;
misdetbound = zeros(1, length(Ss));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xt_s = zeros(3, Ttotal);
phiqt_s = zeros(Ttotal, 1);
thetaqt_s = zeros(Ttotal, 1);
sigmatdb_s = zeros(Ttotal, 1);
betat_s = zeros(Ttotal, 1);
for t = 1 : Ttotal
    xt = zeros(3, 1);
    xt(1:2, 1) = Rmax * 2 * (rand(2, 1)-0.5*ones(2,1));
    while norm(xt) > Rmax || norm(xt) < Rmin || xt(1,1)<0 || abs(atan(xt(2,1)/xt(1,1)))>pi/3
        xt(1:2, 1) = Rmax * 2 * (rand(2, 1)-0.5*ones(2,1));
    end
    xt(3, 1) = 1.5 + (22.5-1.5) * rand(1, 1);
    xt_s(:, t) = xt;
    d2Dt = norm(xt);
    d3Dt = sqrt(d2Dt^2+(-xt(3, 1)+hBS)^2);
    Dt = 4 * (hBS-1) * (xt(3, 1)-1) * fc * 1e9 / c;
    if d2Dt <= Dt
        Plt = 32.4 + 21 * log10(d3Dt) + 20 * log10(fc);
    else
        Plt = 32.4 + 40 * log10(d3Dt) + 20 * log10(fc) - 9.5 * log10(Dt^2 + (hBS - xt(3, 1))^2);
    end
    phiqt_s(t, 1) = rand(1,1) * pi;
    thetaqt_s(t, 1) = rand(1,1) * pi * 0.5;
    sigmatdb_s(t, 1) = -Plt - randn(1, 1) * sigmasfdb;
    betat_s(t, 1) = 10^(0.1*sigmatdb_s(t, 1)) * (randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
for nnids = 1 : length(Ss)
    S = Ss(nnids,1);
    SBS_Basic_init = zeros(Nbsh*Nbsv, S*numberb);
    SBS_Basic_indi_init = zeros(1, S*numberb);
    for p = 1 : Halfnumberb
        for pa = 1 : Halfnumberb
            preal = (p-1) * Halfnumberb + pa;
            SBS_Basic_init(:, (preal-1)*S+1) = SBS(:, preal);
            SBS_Basic_indi_init(1, (preal-1)*S+1) = preal;
        end
    end
    SBS_Basic = SBS_Basic_init;
    SBS_Basic_indi = SBS_Basic_indi_init;
    yqqa = zeros(numberb, Tsetnumber);
    alphaqqa = zeros(numberb, Tsetnumber);
    xtsto = zeros(3, Tsetnumber);
    potenbeamcount = zeros(numberb, 1);
    for t = 1 : Ttotal
        gene_chan_los;
        Propbasic_t;
        if t<=Tsetnumber
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flag = 0;
            for s = 1 : S
                for p = 1 : Halfnumberm
                    for pa = 1 : Halfnumberm
                        aMSx = SMS(:, (p-1) * Halfnumberm + pa);
                        sigmatdb = sigmatdb_s(t, 1);
                        betat = 10^(0.1*sigmatdb) * (randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5);
                        yt_s = abs(betat * aMS'* aMSx)^2 * 10^(MSpower*0.1) /  (10^(0.1*Noisepower)*abs(randn(1,1)*sqrt(2)*0.5 + 1i * randn(1,1)*sqrt(2)*0.5)^2);
                        if yt_s > 10^(0.1*Thetathresdb)
                            flag = 1;
                        else
                        end
                    end
                end
            end
            if flag < 0.5
                misdetection(3, nnids) = misdetection(3, nnids) + 1;
            end
        end
        xxmid = 1;
        phixx = rand(1,1) * pi;
        thetaxx = rand(1,1) * pi * 0.5;
        for nmshcc = 1 : Nmsh
            for nmsvcc = 1 : Nmsv
                aMS((nmshcc-1)*Nmsh+nmsvcc, 1) = 1 / sqrt(Nmsh*Nmsv)...
                    *exp(1i * dvlambda * 2 * pi * ((nmsvcc-1)*sin(phixx)*sin(thetaxx)+(nmshcc-1)*cos(thetaxx)));
            end
        end
        for p = 1 : Halfnumberm
            for pa = 1 : Halfnumberm
                aMSx = SMS(:, (p-1) * Halfnumberm + pa);
                coaMS = abs(aMS' * aMSx)^2;
                u = 10^(0.1*Thetathresdb)*10^(0.1*Noisepower) / (10^(MSpower*0.1)*10^(0.2*sigmatdb)*coaMS);
                fun1=@(y) 1./2.*exp(-0.5.*y).*(1-exp(-0.5.*u.*y));
                xxmid = xxmid * integral(fun1,0,Inf);
            end
        end
        misdetbound(1, nnids) = misdetbound(1, nnids) + xxmid^S;
        disp([nnids, t])
    end
end
misdetection = misdetection / (Ttotal-Tsetnumber);
misdetbound = misdetbound / Ttotal;
%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(Ss, misdetection(4, :), 'k-o','LineWidth',1,'MarkerSize',10)
hold on
semilogy(Ss, misdetbound, 'r-s','LineWidth',1,'MarkerSize',10)
semilogy(Ss, misdetection(3, :), 'b-x','LineWidth',1,'MarkerSize',10)
le = legend( 'Proposed','Bound (16)','Perfect',  'Location', 'northeast');
set(le,'Fontname','Times')
set(gca,'XTick',Ss)
xlabel('Cardinality S','Fontname','Times')
ylabel('PMD','Fontname','Times')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
