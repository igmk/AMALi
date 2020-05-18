function amalifaul(VonDat, BisDat)

% ?ber-Routine zum automatisierten Auswerten der AMALI-Daten
%
% 2/2013
% wertet alle Lidardaten von VonDat bis BisDat aus
%
% VonDat ='170523'; BisDat='170625';     oder '190319' '190408'
% also Format jjmmdd
% Vgl auch Routine superfaul fuer KARL
% 



amalipath % 


campaign = 'AFLUX';          %'ACLOUD'
Speichernamepraefix = 'Version6Flux_';            %'Auswert_BSR_flexausKARL_1s_7_5m';
BSR532soll = 'KARL'; %1.3;
BSRAtFit532start = 1.4;


Starttag=str2num(VonDat(5:6));
Startmonat=str2num(VonDat(3:4));
Startjahr=2000+str2num(VonDat(1:2));
NStart=datenum(Startjahr,Startmonat,Starttag);

Endtag=str2num(BisDat(5:6));
Endmonat=str2num(BisDat(3:4));
Endjahr=2000+str2num(BisDat(1:2));
NEnd=datenum(Endjahr,Endmonat,Endtag);

% UeberlFile= '/home/critter/matlab2/mein/matlab/ueberlappwerte2013_iris4_zpos3';

for j=NStart:NEnd
   tmp = datestr(j); 
   tag = tmp(1:2);
   jahr= tmp(end-1:end);
   if tmp(4:6) =='Jan', monat='01';  kmonat='1'; end
   if tmp(4:6) =='Feb', monat='02';  kmonat='2'; end
   if tmp(4:6) =='Mar', monat='03';  kmonat='3'; end
   if tmp(4:6) =='Apr', monat='04';  kmonat='4'; end
   if tmp(4:6) =='May', monat='05';  kmonat='5'; end
   if tmp(4:6) =='Jun', monat='06';  kmonat='6'; end
   if tmp(4:6) =='Jul', monat='07';  kmonat='7'; end
   if tmp(4:6) =='Aug', monat='08';  kmonat='8'; end
   if tmp(4:6) =='Sep', monat='09';  kmonat='9'; end
   if tmp(4:6) =='Oct', monat='10';  kmonat='a'; end
   if tmp(4:6) =='Nov', monat='11';  kmonat='b'; end
   if tmp(4:6) =='Dec', monat='12';  kmonat='c'; end
   xdatum=[jahr monat tag]; 
   ydatum=[jahr kmonat tag]; 
   
   
   %%% geht es los

   
   disp([upper(mfilename) ': Beginne Standard-Datenauswertung ...']);
   disp(['...vom ', xdatum])
   %ok = amali_eval_ACLOUD(xdatum, campaign, Speichernamepraefix, BSR532soll,BSRAtFit532start);
   ok = amali_eval_Wolke(xdatum, 1, Speichernamepraefix, BSR532soll,BSRAtFit532start);
   
   disp([upper(mfilename) ': Standard-Datenauswertung abgeschlossen.']);
   disp('  ');
    

end




disp([upper(mfilename) ': Das war es dann!!!']);
