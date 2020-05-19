function [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit5( BSRAtFit, FitRange, H, P, PAE, LR, ...
AlRay, BeRa, rueckwaerts)

% Hallo Birte!
  
% Aufruf:
% [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit355arr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
% LR355arr(Sel,i), AlRay355(Sel,i), BeRa355(Sel,i));
% Beta = Beta_total
% Fehler dann:
% BetaerrdP=dBdP.*PAE;
% BetaerrdLR=dBdLR.*15;    bei deltaLR = 15
% BetaerrdR=dBdR.*delta_Randbedg
% CLidar = die sich ergebende Lidarkonstante als Probe
% rueckwaerts: wenn gesetzt wird zum Geraet zurueck integriert


  MFile = [upper(mfilename) ': '];
    
    %%% erste Definition im Falle von Problemen 
    Beta       = P * NaN;
    dBdP     = Beta;
    dBdLR     = Beta;
    dBdR   = Beta; 
    CLidar = Beta;
  wonan=find(isnan(P));
  if length(wonan) == length(P) | length(P) ==0, 
    disp('Signal defekt')
    return
  end
    

FitSel = find(H>FitRange(1) & H<FitRange(2));
beref = BSRAtFit .*mymean(BeRa(FitSel));
X = AlRay-LR.*BeRa;    %Rayleighanteil

if nargin>=8, rueckwaerts = 1; end
if rueckwaerts >0
  
Mz = P.*H.^2 .*exp(-2.*qdrlwvar(H, X));
Nenner= 2.*qdrlwvar(H,LR.*Mz);
Cint = mymean(Mz(FitSel)) ./ beref -mymean(Nenner(FitSel));
Beta = Mz ./ (Nenner + Cint); % Beta-tot

CLidar=  P.*H.^2 ./Beta .*exp(2.*qdrupvar(H,AlRay+LR.*(Beta-BeRa)));
% qdrupvar ist das normale Integral

%dBdR
deltar = 0.01;
beref2 = (BSRAtFit+deltar) .*mymean(BeRa(FitSel));
Cint2 =mymean(Mz(FitSel)) ./ beref2 -mymean(Nenner(FitSel));
Beta2 = Mz ./ (Nenner + Cint2); % Beta-tot
dBdR = (Beta2-Beta)./deltar;  


%dBdLR
deltaLR = 1;
LR2 = LR+deltaLR;
X2 = AlRay-LR2.*BeRa;
Mz2 = P.*H.^2 .*exp(-2.*qdrlwvar(H, X2));
Nenner2= 2.*qdrlwvar(H,LR2.*Mz2);
Cint2 = mymean(Mz2(FitSel)) ./ beref -mymean(Nenner2(FitSel));
Beta2 = Mz2 ./ (Nenner2 + Cint2); % Beta-tot
dBdLR = (Beta2 - Beta) ./ deltaLR;

%dBdP
P2 = P+PAE;
Mz2 = P2.*H.^2 .*exp(-2.*qdrlwvar(H, X));
Nenner2= 2.*qdrlwvar(H,LR.*Mz2);
Cint2 = mymean(Mz2(FitSel)) ./ beref -mymean(Nenner2(FitSel));
Beta2 = Mz2 ./ (Nenner2 + Cint2); % Beta-tot
dBdP = (Beta2-Beta) ./PAE;


else % es soll vorwaerts gerechnet werden

Mz = P.*H.^2 .*exp(2.*qdrupvar(H, X));
Nenner= -2.*qdrupvar(H,LR.*Mz);
Cint = mymean(Mz(FitSel)) ./ beref -mymean(Nenner(FitSel));
Beta = Mz ./ (Nenner + Cint); % Beta-tot

CLidar=  P.*H.^2 ./Beta .*exp(2.*qdrupvar(H,AlRay+LR.*(Beta-BeRa)));


%dBdR
deltar = 0.01;
beref2 = (BSRAtFit+deltar) .*mymean(BeRa(FitSel));
Cint2 =mymean(Mz(FitSel)) ./ beref2 -mymean(Nenner(FitSel));
Beta2 = Mz ./ (Nenner + Cint2); % Beta-tot
dBdR = (Beta2-Beta)./deltar;  


%dBdLR
deltaLR = 1;
LR2 = LR+deltaLR;
X2 = AlRay-LR2.*BeRa;
Mz2 = P.*H.^2 .*exp(2.*qdrupvar(H, X2));
Nenner2= -2.*qdrupvar(H,LR2.*Mz2);
Cint2 = mymean(Mz2(FitSel)) ./ beref -mymean(Nenner2(FitSel));
Beta2 = Mz2 ./ (Nenner2 + Cint2); % Beta-tot
dBdLR = (Beta2 - Beta) ./ deltaLR;

%dBdP
P2 = P+PAE;
Mz2 = P2.*H.^2 .*exp(2.*qdrupvar(H, X));
Nenner2= -2.*qdrupvar(H,LR.*Mz2);
Cint2 = mymean(Mz2(FitSel)) ./ beref -mymean(Nenner2(FitSel));
Beta2 = Mz2 ./ (Nenner2 + Cint2); % Beta-tot
dBdP = (Beta2-Beta) ./PAE;

end

q=1;


  
  
  
  
