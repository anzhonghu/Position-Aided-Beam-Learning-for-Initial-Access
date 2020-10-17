close all;
clear;
LOS_indi = 0;
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
Ncl = 5;
Nray = 20;
ASAdegree = 14.8;
cASAdegree = 22;
Nbsh = 8;
Nbsv = 8;
Nmsh = 4;
Nmsv = 4;
Qb = 6;
Halfnumberb = 2 ^(0.5*Qb);
numberb = 2 ^Qb;
Qm = 4;
Halfnumberm = 2 ^(0.5*Qm);
numberm = 2 ^Qm;
aBS = zeros(Nbsh*Nbsv, 1);
aMS = zeros(Nmsh*Nmsv, 1);
SBS = zeros(Nbsh*Nbsv, 2^Qb);
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
SBS_Basic_init = zeros(Nbsh*Nbsv, S*numberb);
SBS_Basic_indi_init = zeros(1, S*numberb);
for p = 1 : Halfnumberb
    for pa = 1 : Halfnumberb
        preal = (p-1) * Halfnumberb + pa;
        SBS_Basic_init(:, (preal-1)*S+1) = SBS(:, preal);
        SBS_Basic_indi_init(1, (preal-1)*S+1) = preal;
    end
end
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
Npower = [17;20;23;26;29];
access_delay = zeros(2, length(Npower));%2 methods
misdetection = zeros(5, length(Npower));%5 methods
misdetbound = zeros(5, 1);
ahievrate = zeros(5, length(Npower));%5 methods
Tsetnumber = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
for nnids = 1 : length(Npower)
    BSpower = Npower(nnids,1);
    MSpower = Npower(nnids,1);
    succcount = zeros(5,1);
    SBS_Basic = SBS_Basic_init;
    SBS_Basic_indi = SBS_Basic_indi_init;
    yqqa = zeros(numberb, Tsetnumber);
    alphaqqa = zeros(numberb, Tsetnumber);
    xtsto = zeros(3, Tsetnumber);
    potenbeamcount = zeros(numberb, 1);
    for t = 1 : Ttotal
        gene_chan_t;
        if t<=Tsetnumber
        else
            CIin3;
            CIin6;
            Exhaustive;
        end
        Propbasic_t;
        PropSVM_t;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xxmid = 1;
        for p = 1 : Halfnumberm
            for pa = 1 : Halfnumberm
                sigmanx = 0;
                aMSx = SMS(:, (p-1) * Halfnumberm + pa);
                for ncl = 1 : Ncl
                    for nray = 1 : Nray
                        aMS = aMSs(:, (ncl-1)*Nray+nray);
                        coaMS = abs(aMS' * aMSx)^2;
                        sigmanx = sigmanx + sigmanclnraysumfpmd(ncl, 1) * coaMS;
                    end
                end
                if LOSindi > 0.5
                    aMS = aMSs(:,end);
                    coaMS = abs(aMS' * aMSx)^2;
                    sigmany = coaMS * 10^(0.2*sigmatdb);
                    u = 10^(0.1*Thetathresdb)*10^(0.1*Noisepower) /(10^(MSpower*0.1)* (sigmanx * Kr / (Kr+1) + sigmany / (Kr+1)));
                    fun1=@(y) 1./2.*exp(-0.5.*y).*(1-exp(-0.5.*u.*y));
                    xxmid = xxmid * integral(fun1,0,Inf);
                else
                    u = 10^(0.1*Thetathresdb)*10^(0.1*Noisepower) / (10^(MSpower*0.1)*sigmanx);
                    fun1=@(y) 1./2.*exp(-0.5.*y).*(1-exp(-0.5.*u.*y));
                    xxmid = xxmid * integral(fun1,0,Inf);
                end
            end
        end
        misdetbound(nnids, 1) = misdetbound(nnids, 1) + xxmid^S;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([nnids, t])
    end
end
misdetection = misdetection / (Ttotal-Tsetnumber);
misdetbound = misdetbound / Ttotal;




%%%%%%%%%%%%%%%%%%%%%%%%%
plot(Npower, misdetection(1, :), 'k--^','LineWidth',1,'MarkerSize',10)
hold on
plot(Npower, misdetection(2, :), 'k--v','LineWidth',1,'MarkerSize',10)
plot(Npower, misdetection(3, :), 'k--s','LineWidth',1,'MarkerSize',10)
plot(Npower, misdetection(5, :), 'k--x','LineWidth',1,'MarkerSize',10)
plot(Npower, misdetection(4, :), 'k-o','LineWidth',1,'MarkerSize',10)
plot(Npower, misdetbound, 'k-','LineWidth',1,'MarkerSize',10)
le = legend('CI [4]', 'CI [15]', 'Exhaustive [5]', 'SVM [18][19]', 'Proposed','Bound (26)',  'Location', 'northeast');
set(le,'Fontname','Times')
xlim([min(Npower), max(Npower)])
set(gca,'XTick',Npower)
xlabel('BS power (dBm)','Fontname','Times')
ylabel('PMD','Fontname','Times')
grid on%%%%%%%%%%%%%%%%%%%%%%%%
