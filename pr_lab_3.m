clear;
clear all;
f = 3.6 * 10^9;               %częstotliwość
c = 3 * 10^8;                 %prędkość światła
d_max = 2.8;                  %długość pokoju
s_max = 2.6;                  %szerokość pokoju
h_max = 2.6;                  %wysokość pokoju
step1 = 0.005;   
a = -0.5;                       %współczynnik odbicia
h_anten = 1;                  %wysokość anteny
%d = [1:step1:d_max];          %zakres dla zadania 1

d = [0.5 1 1.5 2 2.5 2.8];
absprpo = [];
for i = 1 : 1 : 6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%los
    fTXlosRX(i) = 2 * pi * f * (d(i) ./ c); %LOS
    prpolos(i) = (1 ./ d(i)) .* exp(-1i .* fTXlosRX(1,i));
    absprpolos(1,i) = (abs(prpolos(1,i))) .^ 2;
    absprpo(i,1) = absprpolos(1,i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1-krotne
    dsc(i) = hypot(d(i), s_max);
    dsu(i) = hypot(d(i), (2 * (h_max - h_anten)));
    dpo(i) = hypot (d(i), (2 * h_anten));
    dscp(i) = (2 * d_max) - d(i);
    fTXscRX = 2 * pi * f * (dsc(i) ./ c); %ściana
    fTXsuRX = 2 * pi * f * (dsu(i) ./ c); %sufit
    fTXpoRX = 2 * pi * f * (dpo(i) ./ c); %podłoga
    fTXscpRX = 2 * pi * f * (dscp(i) ./ c); %ściana w linii LOS
    prposc = ((a ./ dsc(i)) .* exp(-1i .* fTXscRX));
    prposu = ((a ./ dsu(i)) .* exp(-1i .* fTXsuRX));
    prpopo = ((a ./ dpo(i)) .* exp(-1i .* fTXpoRX));
    prposcp = ((a ./ dscp(i)) .* exp(-1i .* fTXscpRX));
    absprposc = (abs(prposc))^2;
    absprpo(i,2) = absprposc;
    absprposu = (abs(prposu))^2;
    absprpo(i,3) = absprposu;
    absprpopo = (abs(prpopo))^2;
    absprpo(i,4) = absprpopo;
    absprposcp = (abs(prposcp))^2;
    absprpo(i,5) = absprposcp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2-krotne
    dscsc(i) = hypot(d(i), (2 * s_max)); %ściana + ściana
    dsupo(i) = hypot(d(i), (2 * h_max)); %sufit/podłoga + podłoga/sufit
    dscscp(i) = (2 * d_max) + d(i); % %tył + przód w lini LOS
    fTXscscRX = 2 * pi * f * (dscsc(i) ./ c); 
    fTXsupoRX = 2 * pi * f * (dsupo(i) ./ c); 
    fTXscscpRX = 2 * pi * f * (dscscp(i) ./ c); 
    prposcsc = (((a * a) ./ dscsc(i)) .* exp(-1i .* fTXscscRX));
    prposupo = (((a * a) ./ dsupo(i)) .* exp(-1i .* fTXsupoRX));
    prposcscp = (((a * a) ./ dscscp(i)) .* exp(-1i .* fTXscscpRX));
    absprposcsc = (abs(prposcsc))^2;
    absprpo(i,6) = absprposcsc;
    absprposupo = (abs(prposupo))^2;
    absprpo(i,7) = absprposupo;
    absprposcscp = (abs(prposcscp))^2;
    absprpo(i,8) = absprposcscp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3-krotne
    dscscsc(i) = hypot(d(i), (3 * s_max)); %ściana + ściana + ściana
    dsuposu(i) = hypot(d(i), (2 * (h_max + (h_max - h_anten)))); %sufit + podłoga + sufit
    dposupo(i) = hypot(d(i), (3 * (h_max + h_anten))); %podłoga + sufit + podłoga
    dscscscp(i) = (4 * d_max) - d(i); %tył + przód + tył w lini LOS
    fTXscscscRX = 2 * pi * f * (dscscsc(i) ./ c); 
    fTXsuposuRX = 2 * pi * f * (dsuposu(i) ./ c);
    fTXposupoRX = 2 * pi * f * (dposupo(i) ./ c);
    fTXscscscpRX = 2 * pi * f * (dscscscp(i) ./ c); 
    prposcscsc = (((a * a * a) ./ dscscsc(i)) .* exp(-1i .* fTXscscscRX));
    prpoposupo = (((a * a * a) ./ dsuposu(i)) .* exp(-1i .* fTXsuposuRX));
    prposuposu = (((a * a * a) ./ dposupo(i)) .* exp(-1i .* fTXposupoRX));
    prposcscscp = (((a * a * a) ./ dscscscp(i)) .* exp(-1i .* fTXscscscpRX));
    absprposcscsc = (abs(prposcscsc))^2;
    absprpo(i,9) = absprposcscsc;
    absprpoposupo = (abs(prpoposupo))^2;
    absprpo(i,10) = absprpoposupo;
    absprposuposu = (abs(prposuposu))^2;
    absprpo(i,11) = absprposuposu;
    absprposcscscp = (abs(prposcscscp))^2;
    absprpo(i,12) = absprposcscscp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4-krotne
    dscscscsc(i) = hypot(d(i), (4 * s_max)); %ściana + ściana + ściana + ściana
    dsuposupo(i) = hypot(d(i), (4 * h_max)); %sufit/podłoga + podłoga/sufit + sufit/podłoga + podłoga/sufit
    dscscscscp(i) = (4 * d_max) + d(i); % %tył + przód w lini LOS
    fTXscscscscRX = 2 * pi * f * (dscscscsc(i) ./ c); 
    fTXsuposupoRX = 2 * pi * f * (dsuposupo(i) ./ c); 
    fTXscscscscpRX = 2 * pi * f * (dscscscscp(i) ./ c); 
    prposcscscsc = (((a * a * a * a) ./ dscscscsc(i)) .* exp(-1i .* fTXscscscscRX));
    prposuposupo = (((a * a * a * a) ./ dsuposupo(i)) .* exp(-1i .* fTXsuposupoRX));
    prposcscscscp = (((a * a * a * a) ./ dscscscscp(i)) .* exp(-1i .* fTXscscscscpRX));
    absprposcscscsc = (abs(prposcscscsc))^2;
    absprpo(i,13) = absprposcscscsc;
    absprposuposupo = (abs(prposuposupo))^2;
    absprpo(i,14) = absprposuposupo;
    absprposcscscscp = (abs(prposcscscscp))^2;
    absprpo(i,15) = absprposcscscscp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 15 * length(d);
ds = ([d, dsc, dsu, dpo, dscp, dscsc, dsupo, dscscp, dscscsc, dsuposu, dposupo, dscscscp, dscscscsc, dsuposupo, dscscscscp]);
dss = reshape(ds, 6,15);
Tauds = dss / c;
%Tauds = [Taud Taudsc Taudsu Taudpo Taudscp Taudscsc Taudsupo Taudscscp Taudscscsc Taudsuposu Taudposupo Taudscscscp];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deltaTaudsc = Taudsc - Taud;
% deltaTaudsu = Taudsu - Taud;
% deltaTaudpo = Taudpo - Taud;
% deltaTaudscp = Taudscp - Taud;
% deltaTaudscsc = Taudscsc - Taud;
% deltaTaudsupo = Taudsupo - Taud;
% deltaTaudscscp = Taudscscp - Taud;
% deltaTaudscscsc = Taudscscsc - Taud;
% deltaTaudsuposu = Taudsuposu - Taud;
% deltaTaudposupo = Taudposupo - Taud;
% deltaTaudscscscp = Taudscscscp - Taud;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaTauds = Tauds - Tauds(:,1);
absprpoloss = reshape(absprpolos,[],1);
PrxPr1 = absprpo ./ absprpoloss;
for i = 1 : 1 : 6
    x = deltaTauds(i,:);
    y = PrxPr1(i,:);
    subplot(2,3,i)
    %figure
    stem(x, y)
    set(gca,'yscal','log')
    title(['Profil kanału dla d =  ', num2str(d(i))]);
    xlabel('τ [s]');
    ylabel('Pr/Po');
end
figure
hold on
for i = 1 : 1 : 6
    x = deltaTauds(i,:);
    y = PrxPr1(i,:);
    %subplot(2,3,i)
    %figure
    stem(x, y,'filled')
    set(gca,'yscal','log')
    %title(['Profil kanału dla d =  ', num2str(d(i))]);
    xlabel('τ [s]');
    ylabel('Pr/Po ');
    legend('d = 0.5','d = 1','d = 1.5','d = 2','d = 2.5','d = 2.8')
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : 1 : 6
    AvgTau(i) = (sum(Tauds(:,i) .* PrxPr1(:,i))) / (sum(PrxPr1(:,i)));
    SquareAvgTau(i) = (sum((Tauds(:,i) .^2) .* (PrxPr1(:,i)))) / (sum(PrxPr1(:,i)));
    TauRMS(i) = sqrt((SquareAvgTau(i)) - (AvgTau(i) .^ 2));
end
figure
stem(d,TauRMS)
title('Wykres wartości pierwiastka ze średniokwadratowego rozrzutu opóźnień w funkcji odległości odbiornika od nadajnika');
xlabel('Odległość [m]');
ylabel('τrms [s]');

for i = 1 : 1 : 6
    Bc50(i) = 1 / (5 .* TauRMS(i));
end
figure
stem(d,Bc50)
title('Wykres pasma koherencji');
xlabel('Odległość [m]');
ylabel('Bc50 [Hz]');

for i = 1 : 1 : 6
    Bc90(i) = 1 / (50 .* TauRMS(i));
end
figure
stem(d,Bc90)
title('Wykres pasma koherencji');
xlabel('Odległość [m]');
ylabel('Bc90 [Hz]');