function amalifaul(VonDat, BisDat)

% function to call Amali Processing (2nd step)
%
% 05/2020
% wertet alle Lidardaten von VonDat bis BisDat aus
%
% VonDat ='170523'; BisDat='170625';     oder '190319' '190408'
% Format jjmmdd as string
% 
% 


computer = 'Potsdam';
if strcmp(computer, 'Potsdam')
addpath('/atm_meas/polar_5_6/amali/processing/nadir/cloud/Amali_subroutines');  
end

Speichernamepraefix = 'testtesttest_';            

% reference Value for clean air beneath aircraft. 
% KARL are average values from Ny-Alesund
BSR532soll = 'KARL';%1.2; % %1.3;

%Anfangs(Rate)Wert der unteren Randbedingung vom Klett
BSRAtFit532start = 1.4; 


%formatting strings to fit amali_eval_Wolke
Starttag=str2num(VonDat(5:6));
Startmonat=str2num(VonDat(3:4));
Startjahr=2000+str2num(VonDat(1:2));
NStart=datenum(Startjahr,Startmonat,Starttag);

Endtag=str2num(BisDat(5:6));
Endmonat=str2num(BisDat(3:4));
Endjahr=2000+str2num(BisDat(1:2));
NEnd=datenum(Endjahr,Endmonat,Endtag);


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
   
   %benutzt er ydatum ueberhaupt? 
   %ydatum=[jahr kmonat tag]; 
   
   
   % calling function for 1 File

   disp([upper(mfilename) ': Starting Data processing ...']);
   disp(['...of ', xdatum])
  % ok = amali_eval_Wolkenur532(xdatum, 1, Speichernamepraefix, BSR532soll,BSRAtFit532start,150,155);
   ok = amali_eval_Wolkenur532(xdatum, 1, Speichernamepraefix, BSR532soll,BSRAtFit532start);
   %here aerofak is 1 because we assume aerosol load in NYA and below the
   %aircraft to be the same. if we are aware of deviations from the average
   %in NyA it may me adjusted
   disp([upper(mfilename) '::::::::::: Done ::::::::::' ok]);
   disp('  ');
   disp('  ');
    

end




disp([upper(mfilename) '::::::::::: Finished data processing for given period ::::::::::: ']);
