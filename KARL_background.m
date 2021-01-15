% Jul. 2016 und Jun 15 f?r Hintergrundaerosol MOSAiC da fast keine
% wolkenfreien Daten aus anderen Jahren
%
%
%

clear
% Daten aus:
  LoadDir = '/atm_meas/awipev/lidar/karl/matlab/aer/havg0008/'
% with the name of dattoaer_iterativ_konstantina:
 FileName='_mindata_iterativ.mat';
%FileName='_mindata.mat';

    

SuchStr=[LoadDir '1607' '*' FileName];
%SuchStr=[LoadDir '1706' '*' FileName];
Files = findfile( SuchStr);
if length( Files) == 0,
    disp([MFile 'keine Dateien fuer ' SuchStr ' gefunden!']);
    return;
end
NofFiles = length(Files(:,1));



BSR355arr=[]; 
BSR532arr=[]; 
BSR1064arr=[];
BA355arr=[]; 
BA532arr=[]; 
BA1064arr=[];
BA355Sarr=[]; 
BSR355Sarr=[]; 
BA532Sarr=[]; 
BSR532Sarr=[];
Zeitarr=[];  
VD532arr=[]; 
AA532arr=[];
VD355arr=[]; 
AA355arr=[];
AeroDep532arr=[];
Densarr=[]; 
Temparr=[]; 
O3Densarr=[];
hschritteUVarr=[];
hschritteVISarr=[];
Tau532arr =[];
Tau355arr =[];
Alpha355arr=[];
Alpha532arr=[];

for j=1:NofFiles;
LStr = ['load ' Files(j,:)];
eval(LStr)

BSR355arr=concat(2,BSR355arr,BSR355tot);
BSR532arr=concat(2,BSR532arr,BSR532tot);
BSR1064arr=concat(2,BSR1064arr,BSR1064);
BSR355Sarr=concat(2,BSR355Sarr,BSR355SKlett);
BSR532Sarr=concat(2,BSR532Sarr,BSR532SKlett);
BA355arr=concat(2,BA355arr,BetaAer355tot);
BA532arr=concat(2,BA532arr,BetaAer532tot);
BA1064arr=concat(2,BA1064arr,BetaAer1064);
BA355Sarr=concat(2,BA355Sarr,BetaAer355S);
BA532Sarr=concat(2,BA532Sarr,BetaAer532S);
AA532arr=concat(2,AA532arr,AlphaAer532Raman);
AA355arr=concat(2,AA355arr,AlphaAer355Raman);
VD532arr=concat(2,VD532arr,VD532);
VD355arr=concat(2,VD355arr,VD355);
AeroDep532arr=concat(2,AeroDep532arr,AeroDep532);
%hschritteUVarr=concat(2,hschritteUVarr,hschritteUV);
%hschritteVISarr=concat(2,hschritteVISarr,hschritteVIS);
%Tau355arr=concat(2,Tau355arr,Tau355final);
%Tau532arr=concat(2,Tau532arr,Tau532final);
%Alpha355arr=concat(2,Alpha355arr,Alpha355final);
%Alpha532arr=concat(2,Alpha532arr,Alpha532final);
Zeitarr=concat(2,Zeitarr,DN532);
Densarr=concat(2,Densarr,Density);
Temparr=concat(2,Temparr,Temp);
O3Densarr=concat(2,O3Densarr, O3Density);



end
q=1;

%SuchStr=[LoadDir '1506' '*' FileName];
SuchStr=[LoadDir '1605' '*' FileName];
Files = findfile( SuchStr);
NofFiles = length(Files(:,1));

for j=1:NofFiles;
LStr = ['load ' Files(j,:)];
eval(LStr)

BSR355arr=concat(2,BSR355arr,BSR355tot);
BSR532arr=concat(2,BSR532arr,BSR532tot);
BSR1064arr=concat(2,BSR1064arr,BSR1064);
BSR355Sarr=concat(2,BSR355Sarr,BSR355SKlett);
BSR532Sarr=concat(2,BSR532Sarr,BSR532SKlett);
BA355arr=concat(2,BA355arr,BetaAer355tot);
BA532arr=concat(2,BA532arr,BetaAer532tot);
BA1064arr=concat(2,BA1064arr,BetaAer1064);
BA355Sarr=concat(2,BA355Sarr,BetaAer355S);
BA532Sarr=concat(2,BA532Sarr,BetaAer532S);
AA532arr=concat(2,AA532arr,AlphaAer532Raman);
AA355arr=concat(2,AA355arr,AlphaAer355Raman);
VD532arr=concat(2,VD532arr,VD532);
VD355arr=concat(2,VD355arr,VD355);
AeroDep532arr=concat(2,AeroDep532arr,AeroDep532);
% hschritteUVarr=concat(2,hschritteUVarr,hschritteUV);
% hschritteVISarr=concat(2,hschritteVISarr,hschritteVIS);
% Tau355arr=concat(2,Tau355arr,Tau355final);
% Tau532arr=concat(2,Tau532arr,Tau532final);
% Alpha355arr=concat(2,Alpha355arr,Alpha355final);
% Alpha532arr=concat(2,Alpha532arr,Alpha532final);
Zeitarr=concat(2,Zeitarr,DN532);
Densarr=concat(2,Densarr,Density);
Temparr=concat(2,Temparr,Temp);
O3Densarr=concat(2,O3Densarr, O3Density);



end




crk=NaN*zeros(size(BA532arr));
dimen=size(BA532arr);
wo=find(BSR532arr >1.08 & BSR355arr > 1.04);
crk(wo) =BA355arr(wo) ./ BA532arr(wo);
maxi=zeros(dimen(2),1);
for j=1:dimen(2),
    maxi(j) = max(BSR532arr(1:1400,j)); 
end




% Sommerprofil f?r Birte MOSAiC2020 Flugzeug

dimen2=size(BSR532arr,2);
niederzahl=40;
gutposi=[];
ber=5:600;
for j=1:dimen2,
    z1=find(BSR532arr(ber,j)< 1.0);
    z1l=length(z1);
    z2=max(BSR532arr(ber,j));
    z3=max(VD532arr(ber,j));
    if z1l < niederzahl & z2 <2.5 & z3<0.02 , gutposi=concat(1,gutposi,j); end
end

BSR532mean=(mean(BSR532arr(:,gutposi)'))';
BSR532Smean=(mean(BSR532Sarr(:,gutposi)'))';
BSR355mean=(mean(BSR355arr(:,gutposi)'))';
BSR532median=H.*0;
BSR532Smedian=H.*0;
BSR355median=H.*0;
for j=1:length(H);
    BSR532median(j)=median(BSR532arr(j,gutposi));
    BSR532Smedian(j)=median(BSR532Sarr(j,gutposi));
    BSR355median(j)=median(BSR355arr(j,gutposi));
end


%save /atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2020/aerosol_background_karl/KARLaverageBSR532ausSommer H BSR532median ...
%    BSR532Smedian BSR355median BSR532mean BSR532Smean BSR355mean 