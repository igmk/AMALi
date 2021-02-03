function ok = amali_eval_Wolkenur532(DatStr, aerofak, Speichernamepraefix, BSR532soll, BSRAtFit532start, Von, Bis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AMALi Auswertung
%
% fuer senkrecht nach unten und Iteration in Wolken
%
% vlg altes amali_eval_ACLOUD hier eine Zeititeration f?r alle Kanaele
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DatStr='170523'
% aerofak=1; %Faktor um NyA Aerosol R?ckstreuung an den aktuellen AMALi
% Flug zu skalieren
% Speichernamepraefix = 'Tollertest'
% BSR532soll = 1.3;  angestrebter Referenzwert von BSR532 unter klaren
% Bedingungen   wenn BSR532soll = 'KARL' dann liest er den KARL-NyA Wert
% ein
% BSRAtFit532start   Anfangs(Rate)Wert der Randbedingung vom Klett
% Von / Bis: optional Die Datensatznummern in denen gerechnet werden soll
% wird Von / Bis nicht gesetzt rechnet er alles
% ok = amali_eval_Wolke('170523', 1, 'Test1', 1.3, 1.4, 2000, 2100)
% ok = amali_eval_Wolke('170523', 1, 'Test1', 'KARL', 1.4)




% Klett wird nur fuer 532 duchgefuerht. Hintergrundkorrektur und korrektur
% wegen offset wird fuer alle kanaele durchgefuehrt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ok=-1;

%
disp("Birtes Version")


% Einstellungen fuer Kampagnen
% depending on the date it looks up campaign and input data
% has to be adjusted if run in Cologne

%Angstroem exponent relevant for conversion of input between wavelength
% not relevant if only 532 is calculated
computer = 'Potsdam';
if strcmp(computer, 'Potsdam')
    if strcmp(DatStr(1:2), '17')
        campaign = 'ACLOUD';
        Angstroem=1.4; %#ok<NASGU>
        speicherdir = '/atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2017/';
        amalidatendir='/lidar4/lidar/amali2011/data/mat/2017/'
    elseif strcmp(DatStr(1:2), '19')
        campaign = 'AFLUX';
        Angstroem=1.2; %#ok<NASGU>
        speicherdir = '/atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2019/';
        amalidatendir='/atm_meas/polar_5_6/amali/data/mat/2019/';
        
    elseif strcmp(DatStr(1:2), '20')
        campaign = 'MOSAICACA';
        Angstroem=1.2; %#ok<NASGU>
        speicherdir = '/atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2020/';
        amalidatendir='/atm_meas/polar_5_6/amali/data/mat/2020/';
    else
        disp('Campaign unknown')
        return
    end
end
%only needed if input is converted between wavelengths
%WvlExpo = 4-Angstroem;

% wenn nur ein teil der daten gerechnent wird, sonst kann von:bis auch leer sein - dann nimmt er alles
if nargin > 5
    teilzeitjob = 1;
else
    teilzeitjob = 0;
end

% Definitionen ------------------------------------------------------------
Wvl532 =5.3207e-7;
Wvl355 =3.5471e-7;


% LR = lidar ratio = attenuation/backscatter
% the smaller the droplets the larger the LR (fig. 3 in O'Connor 2004)
LR532aerosol=23; %(clean marine air KIM 2018)
LR532wolke=15; %(Hogan 2003 LR for Liquid with effective radii 2-5mym) % YORKS 2011 8-60 eher 11 for low stratiform water clouds, 17 for higher water clouds)
LR532wolke2=23 % as clean air, not sensitive to wrong assumption in cloud mask, in the range of realistic cloud LR
LR532wolke3= 27 %(LAMPERT 2009 for thin ice cloud in Arctic)

LRerr = 2; % angenommener Fehler im LR %fuer das ausrechnen von dbeta


Wolkenschwelle532=5;% wenn das bsr (backscatter ratio = gemessener backscatter / ausgerechneter molekularer Backscatter) groesser als 5 ist nehmen wir an, dass es sich um eine Wolke handelt und nutzen das LR fuer wolke


%tzlim = 5;
hzlim=5;  %fuer die while schleifen, limit f?r Hoch- Tiefzahl bei Iterationenen, so oft versucht er es nochmal wenn LR auserhalb des erlaubten Berreiches liegt

% BSR at LBC (lower boundary condition)
% LBC between 7.5 und 107.5 m ueber ground unter der annahme, das Amali senkrecht guckt
% BSRAtFit532start = 1.4;      % besser als Inputparameter
% nur relevant wenn auch 355 gerechnet werden soll:
% BSRAtFit355start = 1+(BSRAtFit532start-1) ./ 1.5.^WvlExpo; % rechnen 355er Randbedingung aus der von 532

%minimum BSR at UBC (upper boundary condition)
%UBC zwischen ende des ueberlaps und der ersten Wolke
BSR532mintrust  = 1.2; %ist das jetzt das minimum oder das maximum, dass wir akzeptieren?

%BSR355mintrust  = 1 + (BSR532mintrust -1) ./ 1.5.^WvlExpo;

% BSR532soll = 1.3;  % ist Inputparameter
%nimmt er aus dem Karl, es soll da rauskommen das die Aerosolbelastung wie
%in NyA ist

% nur f?er den Notfall (NaN im KARL oder keine guteaeroposi) fuer UBC:
BSR532sollnotfall = 1.3;
%BSR355sollnotfall = 1 + (BSR532sollnotfall-1) ./ 1.5.^WvlExpo;

%300 m D
UeberlappEnde=300; % [m]
%pretrigger = Aufnahme des Signals bevor der Laser abgeschossen wird
% noetig um Hintergrund und Hintergrundrauschen zu bestimmen
pretrigrange=1:400; %in dem bereich bestimmen wir das Rauschen
pretriggerbins=405; %hier faengt das eigentliche Signal an

Schwelle = 1e-8; % Signal (P)< Schwelle nur noch Rauschen

UeberlappEndeposi = round(UeberlappEnde ./ 7.5); %von m in bin umrechnen

diffisoll = 0.05;  % so gut wollen wir BSR bestimmen
%FitRange=[2600, 2700]; %Standardwahl, wenn es nichts Besseres gibt

Hcalcrange =3500;    %bis 3500m entfernung vom Geraet rechnet er.
%muss hoch gesetzt werden wenn das Flugzeug ueber 3500 m fliegt



% ------------------------------------------------------------------------
% Rayleighstreuung aus median Radiosonde NyA, laedt meadian KARL
if strcmp(computer, 'Potsdam')
    if strcmp(campaign, 'MOSAICACA')
        ptuinfile='/atm_meas/awipev/lidar/karl/matlab/ptu/2009.mat';
        load(ptuinfile)
        %     ozoinfile='/atm_meas/awipev/lidar/karl/matlab/ozo/2009.mat';
        %     load(ozoinfile)
        load /atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2020/aerosol_background_karl/KARLaverageBSR532ausSommer
        KarlH = H; % H wird sp?ter der AMALi Range-vektor
        %BSR532Karlmean= BSR532mean;
        BSR532Karlmedian = BSR532median.*aerofak;
        %BSR532SKarlmean= BSR532Smean; BSR532SKarlmedian = BSR532Smedian;
        %BSR355Karlmean= BSR355mean; BSR355Karlmedian = BSR355median;
        load /atm_meas/polar_5_6/flight_data/gps/MOSAICACAGPS.txt -ascii
        flugzeit=MOSAICACAGPS(:,1);
        flughoehe=MOSAICACAGPS(:,2);
        wo=find(flughoehe > 10000);
        flughoehe(wo)=NaN;
    elseif strcmp(campaign, 'ACLOUD')
        ptuinfile='/atm_meas/awipev/lidar/karl/matlab/ptu/1706.mat';
        load(ptuinfile)
        %     ozoinfile='/atm_meas/awipev/lidar/karl/matlab/ozo/1706.mat';
        %     load(ozoinfile)
        load /atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2017/aerosol_background_karl/KARLaverageBSR532ausACLOUD
        KarlH = H; % H wird sp?ter der AMALi Range-vektor
        %BSR532Karlmean= BSR532mean;
        BSR532Karlmedian = BSR532median.*aerofak;
        %BSR532SKarlmean= BSR532Smean; BSR532SKarlmedian = BSR532Smedian;
        %BSR355Karlmean= BSR355mean; BSR355Karlmedian = BSR355median;
        load /atm_meas/polar_5_6/flight_data/gps/ACLOUDGPS.txt -ascii
        flugzeit=ACLOUDGPS(:,1);
        flughoehe=ACLOUDGPS(:,2);
        wo=find(flughoehe > 10000);
        flughoehe(wo)=NaN;
    elseif strcmp(campaign, 'AFLUX')
        ptuinfile='/atm_meas/awipev/lidar/karl/matlab/ptu/1904.mat';
        load(ptuinfile)
        %was macht er hier eigentklich im Maerz? Ist das tragisch?
        %--- ggf. Christoh fragen
        %     ozoinfile='/atm_meas/awipev/lidar/karl/matlab/ozo/1904.mat';
        %     load(ozoinfile)
        load /atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2019/aerosol_background_karl/KARLaverageBSR532ausAFLUX
        KarlH = H; % H wird sp?ter der AMALi Range-vektor
        %BSR532Karlmean= BSR532mean;
        BSR532Karlmedian = BSR532median.*aerofak;
        %BSR532SKarlmean= BSR532Smean; BSR532SKarlmedian = BSR532Smedian;
        %BSR355Karlmean= BSR355mean; BSR355Karlmedian = BSR355median;
        load /atm_meas/polar_5_6/flight_data/gps/AFLUXGPS.txt -ascii
        flugzeit=AFLUXGPS(:,1);
        flughoehe=AFLUXGPS(:,2);
        wo=find(flughoehe > 10000);
        flughoehe(wo)=NaN; %#ok<*FNDSB>
    end
end
%Ozon wird ignoriert da unter 4000m Absorption durch Ozon nicht relevant
%meanO3profile=(mymean(OZOO3Density'))';

%calculate Air density
Density = density(PTUPressure, PTUPressure./100, PTUTemperature,ones(size(PTUTemperature)).*2); %#ok<*NODEF>
%Extinktion nur Streuanteil (99% ausmacht)
%alpha ist gleich fuer beide polarisationen
PTUAlRay532 = Density .* raytotwq( Wvl532, PTUTemperature, Density);
PTUAlRay355 = Density .* raytotwq( Wvl355, PTUTemperature, Density);

% Rayleigh R?ckstreuung (die schwach polaris.- und temp- abhaengig ist
% Streuung =  Luftdichte (Density) mal Wirkungsquerschnitt
PTUBeRa532=Density.*raybckwq (Wvl532,'p','p', PTUTemperature, Density);
PTUBeRa532S=Density.*raybckwq (Wvl532,'p','s', PTUTemperature, Density);
PTUBeRa355=Density.*raybckwq (Wvl355,'p','u', PTUTemperature, Density);
%Warum ist Beta eigentlich polarisationsabhaengig??? wir gucken doch auf
%zufaellig ausgerichtete Profile


% Korrektur durch molek. Absorption, die praktisch nur durch O3
% vor allem wichig in Strato - hier vermutlich egal
% ptudimen=size(PTUAlRay355);
% for j=1:ptudimen(2)
%     ARayab532=meanO3profile.*o3abswq(Wvl532,PTUTemperature(:,j));
%     ARayab355=meanO3profile.*o3abswq(Wvl355,PTUTemperature(:,j));
%     PTUAlRay532(:,j)=PTUAlRay532(:,j)+ARayab532;
%     PTUAlRay355(:,j)=PTUAlRay355(:,j)+ARayab355;
% end


% gemittelte Profile aus Ny-Alesund
% ist es schlau hier ein Mittel zu nehmen? Vielleicht beesser daten des
% jeweiligen Tags? Ergebnisse schwanken um etwa 7%
%rayleigh profile werden auf allen kanaelen mitgerechnet
% noch Ny_Alesund H?henvektor
NAlRay532=(mymean(PTUAlRay532'))';
NAlRay355=(mymean(PTUAlRay355'))';
NBeRa532=(mymean(PTUBeRa532'))';
NBeRa355=(mymean(PTUBeRa355'))';
NBeRa532S=(mymean(PTUBeRa532S'))';

NH = PTUHeight(:,1);
NDensity=(mymean(Density'))';

% In den ersten bins aus Ny-Alesund kann NaN stehen, dies interpolieren wir
% weg
tmp=real(log(NAlRay532)); %gerade kommt raus, weil Rayleigh profil quasi exponentiell (folgt Luftdichte)
wo=find(isnan(tmp(1:20)));
if ~isempty(wo)
    woend=wo(end);
    st=mymean(diff(tmp(woend+1:woend+11)));
    for jj=woend:-1:1
        tmp(jj)=tmp(woend+1)-(woend+1-jj).*st;
    end
end
NAlRay532 = exp(tmp);

tmp=real(log(NBeRa532));
wo=find(isnan(tmp(1:20)));
if ~isempty(wo)
    woend=wo(end);
    st=mymean(diff(tmp(woend+1:woend+11)));
    for jj=woend:-1:1
        tmp(jj)=tmp(woend+1)-(woend+1-jj).*st;
    end
end
NBeRa532 = exp(tmp);


%-------------------------------------------------------------------------


% Einlesen der Amali-Lidardaten


% Fuer alte Datenstruktur: datum=[DatStr '/'];
%aus unerfindlichen gruenden haben die neuen matlabfiles keine nullen vor dem Monat

if strcmp(campaign, 'ACLOUD')
    datum=[DatStr '/'] %#ok<*NOPRT>
else
    datum=[DatStr(1:2), DatStr(4:6)];
end
suchfile=[amalidatendir datum '*.mat'];

files=findfile(suchfile)
NofFiles = length(files(:,1));
if NofFiles >1
    disp('mehr als einen Datensatz gefunden, alles wird durchiteriert')
end


if isempty(files)
    ok=-1;
    disp(['.... ... keine Daten fuer: ' datum])
    return
else
    %
    ok = NofFiles;
    for kk=1:NofFiles %Schleife ?ber alle Datens?tze
        disp(['working on datasset ' num2str(kk) ' of ' num2str(ok)])
        load(files(kk,:))
        %falls bei den Dateinamen was komisch ist
        wostrich=strfind(files(kk,:),'/');
        if length(wostrich>1), wostrich=wostrich(end); end %#ok<ISMT>
        womat= strfind(files(kk,:), '.mat');
        filename = files(kk,wostrich+1:womat-1);
        
        
        datasize=size(alldata);
        entries=datasize(2);
        bins=datasize(3);
        data532=(reshape(alldata(1,:,:),entries, bins));   % Format umgedreht jetzt wie KARL
        
        %andere Kanaele werden auch eingelesen um Hitnergrundkorrektur
        %durchzufuehren. und am Ende das Signal im selben format wie fuer
        %532 wegzuschreiben
        data532s=(reshape(alldata(3,:,:),entries, bins));  % size =[H, Zeiten]
        data355=(reshape(alldata(5,:,:),entries, bins));
        
        %counting kanaele (nutzen wir nicht)
        %data532c=(reshape(alldata(2,:,:),entries, bins));
        %data532sc=(reshape(alldata(4,:,:),entries, bins));
        %data355c=(reshape(alldata(6,:,:),entries, bins));
        
        
        matlabzeit = (allinfo(1,:)/86400) + datenum(1904,1,1);
        
        H=(0:7.5:3500)';    %dH=7.5;    lim=50; %(Limit Signal hat Wolken) %??? keine Ahnung mehr wafuer er das genau benutzt
        
        LH=length(H);
        P532A=zeros(entries,LH);
        P532SA=zeros(entries,LH);
        P355A=zeros(entries,LH);
        
        
        % Im ersten Datensatz (k=1) wandeln wir die Karl H?hen und BSRWerte auf das
        % AMALi Gitter um
        if kk==1
            if BSR532soll=='KARL'
                [KarlHzumBoden,indexx]=sort(KarlH,'descend');
                BSR532KarlmedianzumBoden = BSR532Karlmedian(indexx);
                % jetzt m?ssen wir das auf das H?hengitter des Amali interpolieren
                BSR532Karlmedianvgl = interp1(KarlHzumBoden, BSR532KarlmedianzumBoden, H);
                
            else
                
                BSR532Karlmedianvgl = repelem(BSR532soll,length(H))';
            end
        end
        
        
        %Hintergrundkorrektur
        P532Abackground=zeros(entries,1);
        P532SAbackground=zeros(entries,1);
        P355Abackground=zeros(entries,1);
        
        
        %spatial offsets - personal communication Roland Neuber, Checked so
        %that cloud tops align.
        %         vs532p=-1;
        %         vs532s=-1;
        %         vs355=+1;
        %nicht loeschen, wichtige info!!!! die kaneale passen nicht genau
        %aufeinander
        %wenn jetzt nicht alle Kanaele prozessiert werde muss das hinterher
        %drauf gerechnet werden.
        %ggf. Kanaele trotzdem mit einlesen und P verschieben und ausschreiben??
        %Wie wissen wir das diese Info im Absoluten Raum korrekt ist ???
        %Vergleich Bodensignal???
        % Porblem besteht nur fuer counting signal ----
        
        %auch nur counting signal -- also hier ertmal egal
        %wo = find(data532c > 3.2e4);
        %data532creserve=data532c;
        %data532c(wo)=NaN;
        %wo = find(data532sc > 3.0e4);
        %data532screserve=data532sc;
        %data532sc(wo)=NaN;
        %wo = find(data355c > 3.2e4);
        %data355creserve=data355c;
        %data355c(wo)=NaN;
        
        
        % Hintergrund wird berechnet und vom Signal abgezogen
        for j=1:entries
            P532Abackground(j) = mymean(data532(j,pretrigrange));
            P532A(j,:)=data532(j,pretriggerbins+1:pretriggerbins+LH)- P532Abackground(j);
            P532SAbackground(j) = mymean(data532s(j,pretrigrange));
            P532SA(j,:)=data532s(j,pretriggerbins+1:pretriggerbins+LH)-P532SAbackground(j);
            P355Abackground(j) = mymean(data355(j,pretrigrange));
            P355A(j,:)=data355(j,pretriggerbins+1:pretriggerbins+LH)- P355Abackground(j);
            %P532C(j,:)=data532c(j,vs532p+pretriggerbins+1:vs532p+pretriggerbins+LH)-mymean(data532c(j,pretrigrange));
            %P532SC(j,:)=data532sc(j,vs532s+pretriggerbins+1:vs532s+pretriggerbins+LH)-mymean(data532sc(j,pretrigrange));
            %P355C(j,:)=data355c(j,vs355+pretriggerbins+1:vs355+pretriggerbins+LH)-mymean(data355c(j,pretrigrange));
            
        end
        
        %rechnen mit Lidarsigal P Klett --> transponiert, damit das format
        %ist wie bei Karl
        P532Klett =P532A';
        P532SKlett =P532SA';
        P355Klett =P355A';
        
        
        %unnoetig, kann im naechsten schritt berechnen und kalibriert werden
        %         P532fuerdepol= P532A';
        %         P532Sfuerdepol= P532SA';
        %         VolDep532 = P532Sfuerdepol ./ P532fuerdepol;
        
        %Wenn erte negativ oder null werden werden sie durch sehr kleinen
        %Wert (Schwelle) ersetzt
        wo=find(P532Klett < Schwelle); P532Klett(wo)=Schwelle;
        wo=find(P532SKlett < Schwelle); P532SKlett(wo)=Schwelle;
        wo=find(P355Klett < Schwelle); P355Klett(wo)=Schwelle;
        
        P532Klettnoise=zeros(size(P532Klett));
        backgroundnoise532=zeros(entries,1);
        P532SKlettnoise=zeros(size(P532SKlett));
        backgroundnoise532S=zeros(entries,1);
        P355Klettnoise=zeros(size(P355Klett));
        backgroundnoise355=zeros(entries,1);
        
        % warum ausgerechnet 1.5% mal die Std ist mir raetselhaft
        for j=1:entries
            backgroundnoise532(j)=std(data532(j,pretrigrange)).*1.5;
            P532Klettnoise(:,j)=real(sqrt(P532Klett(:,j)))+backgroundnoise532(j);
            backgroundnoise532S(j)=std(data532s(j,pretrigrange)).*1.5;
            P532SKlettnoise(:,j)=real(sqrt(P532SKlett(:,j)))+backgroundnoise532S(j);
            backgroundnoise355(j)=std(data355(j,pretrigrange)).*1.5;
            P355Klettnoise(:,j)=real(sqrt(P355Klett(:,j)))+backgroundnoise355(j);
        end
        SNR532=P532Klett./P532Klettnoise; %#ok<NASGU>
        SNR532S=P532SKlett./P532SKlettnoise; %#ok<NASGU>
        SNR355=P355Klett./P355Klettnoise; %#ok<NASGU>
        
        
        % Rayleighgr?ssen + Density auf Amaligitter interpolieren
        NAlRay532aufl=interp1(NH,NAlRay532,H);
        AlRay532=NAlRay532aufl(end:-1:1);
        NAlRay355aufl=interp1(NH,NAlRay355,H);
        AlRay355=NAlRay355aufl(end:-1:1); %#ok<NASGU>
        NBeRa532aufl=interp1(NH,NBeRa532,H);
        BeRa532=NBeRa532aufl(end:-1:1);
        NBeRa532Saufl=interp1(NH,NBeRa532S,H);
        BeRa532S=NBeRa532Saufl(end:-1:1); %#ok<NASGU>
        NBeRa355aufl=interp1(NH,NBeRa355,H);
        BeRa355=NBeRa355aufl(end:-1:1); %#ok<NASGU>
        NDensityaufl=interp1(NH,NDensity,H);
        Density=NDensityaufl(end:-1:1);
        
        
        % Nun der schlaue Klett
        % wir brauchen f?r jeden Kanal:
        % BSRAtFitarr = das array der BSR-Werte f?r den Fit
        % FitRange = die H?he von, bis in der im Durchschnitt BSR=BSRAtFitarr sein
        % soll z.B. [2000, 2400] % Upper boundary condition???
        % H: Hoehenvektor
        % P: das Lidarprofil
        % Perr "error" Fehler des Lidarprofils (=Rauschen/Noise)
        % LRarr= array des Lidarverh?ltnisses, gleich gro? wie P
        % ALRay532 / BeRa532 die Rayleigh-Extinktions und R?ckstreuwerte
        % klettinv_ableit4 berechnet:
        % 1) Beta = Gesamtr?ckstreuung (aer+Rayleigh)
        % 2) dBdR die partielle Ableitung d Beta_aer / d Randbedingung (oben oder
        % unten???
        % 3) dBdLR= die part. Ableit. d Beta_aer / d Lidarverh?eltnis
        % 4) dBdP: die part. Abl. d Beta_aer / d Perr(rauschen)
        % 5) Vorschlag CLidar: die Lidarkonstane aus L?sung Beta, LR und elast.
        % Lidargleichung
        
        % Definitionen:
        % f?r alle Farben:
        FitRangearr=zeros(2,entries);
        Wolkenmaske=zeros(size(P532Klett)); %
        dimen=size(P532Klett);
        
        % f?r 532P
        LR532arr=ones([size(P532Klett),3]).*LR532wolke;    % es geht ja um Wolken
        LR532arrerr=ones([size(P532Klett),3]).*LRerr;
        BSRAtFit532arr=ones(entries,1).* BSRAtFit532start;
        Btemp532=zeros(LH,1); Btemp532_2=zeros(LH,1);    Btemp532err=Btemp532;
        Betaaer532=zeros(LH,1); Betaaer532err=Btemp532;  BetaAer532_2=zeros(LH,1);
        BetaAer532Klett=zeros(size(P532Klett)); BetaAer532Kletterr=zeros(size(P532Klett));
        %BetaAer532Klettfest=zeros(size(P532Klett)); BetaAer532Klettfesterr=zeros(size(P532Klett));
        dBeta532dLR=zeros(size(P532Klett));
        dBeta532dR=zeros(size(P532Klett));
        dBeta532dP=zeros(size(P532Klett));
        BSR532Klett=zeros([size(P532Klett),3]); BSR532Kletterr=zeros(size(P532Klett));
        %BSR532Klettfest=zeros(size(P532Klett)); BSR532Klettfesterr=zeros(size(P532Klett));
        C532Lidar=zeros([size(P532Klett),3]);
        AlphaAer532= zeros([size(P532Klett),3]);
        AlphaAer532err= zeros([size(P532Klett),3]);
        attenu532P= zeros([size(P532Klett),3]);       attenu532Perr= zeros(size(P532Klett));
        Abbruch532 = zeros(dimen(2),3);
        BSR532sollarr = zeros(dimen(2),3);
        % Abbruch532 gibt Abbruchkriterium an:
        % 0:hat gar nicht LR angepa?t
        % 1: normale Konvergenz
        % 2: LR = 5 (unteres Limit)
        % 3: LR = LRobergr (oberes Limit)
        % 4: hat sich totiteriert ohne Konvergenz
        % 5: guteaeroposi = 0;
        %
        
        
        Ptheo532=2.*BeRa532./H.^2.*exp(-2.*qdrupvar(H,2.*AlRay532)); %was passiert hier? Wo kommt die 2 her? ansonsten ist das die Lidargleichung fuer reine rayleigh streuung
        
        % es hat etwas zu tun das wir das signal erstmal auf die richtige
        % groessenordnug skalieren, aber was genau - ???
        
        
        
        
        
        
        if teilzeitjob==1
            start = Von; ende = Bis;
        else
            start = 1; ende = entries;
        end
        teiler=floor(ende./20);
        
        
        %schleife ?ber Lidar ratios
        for LR=1:3
            if LR ==1
                LR532arr(:,:,LR)=ones(size(P532Klett)).*LR532wolke;
            elseif LR ==2
                LR532arr(:,:,LR)=ones(size(P532Klett)).*LR532wolke2;
            elseif LR ==3
                LR532arr(:,:,LR)=ones(size(P532Klett)).*LR532wolke3;
            else
                'number of LR not correct'
            end
            
            % Hier geht es los: die gro?e Schleife ueber alle Zeitschritte
            %
            %
            
            for j=start:ende   %1:entries    %2000:5000  start:ende
                
                if (mod(j,teiler)) ==0, disp(['Profile: ' num2str(j) ' von ' num2str(ende)]), end
                if j ==5000
                    q=1;%#ok<NASGU> %???
                end
                
                
                % Definitionen innerhalb eines Zeitschrittes
                Psoll532=Ptheo532 ./Ptheo532(UeberlappEndeposi).*P532Klett(UeberlappEndeposi,j); % ???  ich hab keinen Schimmer was hier passiert
                %Psoll532S=Ptheo532S ./Ptheo532S(UeberlappEndeposi).*P532SKlett(UeberlappEndeposi,j);
                %Psoll355=Ptheo355 ./Ptheo355(UeberlappEndeposi).*P355Klett(UeberlappEndeposi,j);
                ichmerkmirkomischepositionen = 0; %#ok<NASGU>
                
                
                
                % also da wo
                % "select" die guten H?henbins, in denen gerechnet werden soll.
                %
                Sel532P = connrnge( H >= 0 & H <= Hcalcrange & ...
                    P532Klett(:,j) > 0 & ...
                    Density(:,1) > 0, 1);
                if length(Sel532P) > 1
                    Sel532P = (Sel532P(1):Sel532P(2));
                else
                    Sel532P = [];
                end
                
                %             % "select" die guten H?henbins, in denen gerechnet werden soll.
                %             %
                %             Sel532S = connrnge( H >= 0 & H <= Hcalcrange & ...
                %                 P532SKlett(:,j) > 0 & ...
                %                 Density(:,1) > 0, 1);
                %             if length(Sel532S) > 1
                %                 Sel532S = (Sel532S(1):Sel532S(2));
                %             else
                %                 Sel532S = [];
                %             end
                %
                %
                %             % "select" die guten H?henbins, in denen gerechnet werden soll.
                %             %
                %             Sel355P = connrnge( H >= 0 & H <= Hcalcrange & ...
                %                 P355Klett(:,j) > 0 & ...
                %                 Density(:,1) > 0, 1);
                %             if length(Sel355P) > 1
                %                 Sel355P = (Sel355P(1):Sel355P(2));
                %             else
                %                 Sel355P = [];
                %             end
                %
                
                
                
                
                
                
                %set range gates for lower boundary condition
                [~,posi]=min(abs(matlabzeit(j)-flugzeit));
                hwo2=flughoehe(posi)-7.5;
                hwo1=flughoehe(posi)-107.5;
                FitRangearr(1,j)=hwo1; FitRangearr(2,j)=hwo2;
                
                
                
                
                % wir rechnen je ein erstes Mal, nur um zu sehen, wo Wolken sind
                % wenn FitRangearr(:,j) ungefaehr gleich Hcalcrange, dann ist das
                % Lidarsignal durch die Wolke durchgekommen, also ist die Wolke d?nn und
                % damit hat BSRAtFit einen Einflu?. Dann k?nnte Beta zu klein werden
                % wir iterieren BSRAtFit532arr(j) bis Beta nicht mehr offensichtlich zu
                % klein ist und benutzen median(btemp) > 1.05 als Abbruchkriterium
                % Sinn: eine Wolkenmaske zu erstellen
                
                
                % 532P
                condi=1; iter=0; itmax=600;
                hz532=0;
                
                %finde den Berreich als 'clear '
                %F?r die erste iteration wird der der Bereich der UBC einfach
                %festgenagelt und die LBC solange erh?ht bis BSR532mintrust
                %eerreicht wird. dient einer groben scalierung um dann im
                %n?chsten Schritt die Wolkenmaske abzuleiten
                [VoI,BiI]= minintsuche(abs(Psoll532(Sel532P(40:200))-P532Klett(Sel532P(40:200),j)),40);
                clearint=VoI:BiI;
                while condi
                    iter = iter+1;
                    if iter > itmax, condi = 0; end
                    
                    BSRAtFiterr = BSRAtFit532arr(j) ./ 5; %warum wird hier durch 5 geteilt ????  irgendwas fuer fehlerabschaetzung
                    
                    [Beta, dBdR532, dBdLR532, dBdP532, CLidar] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Klettnoise(Sel532P,j), ...
                        LR532arr(Sel532P,j,LR), AlRay532(Sel532P,1), BeRa532(Sel532P,1)); %#ok<ASGLU>
                    Betaaer532(Sel532P)=Beta-BeRa532(Sel532P,1);
                    Betaaer532err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Klettnoise(Sel532P,j));
                    Btemp532(Sel532P)=Beta./BeRa532(Sel532P,1); Btemp532err(Sel532P)=Betaaer532err(Sel532P,j)./BeRa532(Sel532P,1);
                    tmp=mymean(Btemp532(Sel532P(clearint)));
                    
                    [Beta2, ~, ~, ~, ~] = klettinv_ableit4( BSRAtFit532arr(j)+1, FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Klettnoise(Sel532P,j), ...
                        LR532arr(Sel532P,j,LR), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                    BetaAer532_2(Sel532P)=Beta2-BeRa532(Sel532P,1);
                    %BetaAer532_2err(Sel532P,j)=abs(dBdR5322.*BSRAtFiterr)+abs(dBdLR5322.*LR532arrerr(Sel532P,j))+abs(dBdP5322.*P532Kletterr(Sel532P,j));
                    Btemp532_2(Sel532P)=Beta2./BeRa532(Sel532P,1); %Btemp532_2err(Sel532P)=BetaAer532_2err(Sel532P,j)./BeRa532(Sel532P,1);
                    tmp2=mymean(Btemp532_2(Sel532P(clearint)));
                    
                    deltasoll = BSR532mintrust - tmp;
                    deltab=(tmp2-tmp);
                    
                    %Zum plotten der zwischenschritte
                    %                 BetaKlett1(:,iter,j)= Beta;
                    %                 LBCKlett1(:,j)=BSRAtFit532arr(j);
                    %                 UBCKlett1(:,iter,j)=  Btemp532;
                    %                 CKlett1(:,iter,j)=CLidar;
                    %                 clearKlett1(:,iter,j)=clearint;
                    %                 HKlett1(:,j)=H(Sel532P);
                    
                    if abs(deltasoll) > 0.01 %Ziel noch nicht erreicht
                        %w=(BSR532sollarr(j) - tmp) ./ deltab;
                        BSRAtFit532arr(j) = BSRAtFit532arr(j) +1./deltab.*deltasoll; %Passe LBC an
                        if BSRAtFit532arr(j) > 1e6, hz532=hz532+1; BSRAtFit532arr(j) =1e6; end %Wenn LBC super unrealistisch hoch setze es wieder auf einen Wert der halbwegs realistisch ist (Daamit es nicht ausversehen irgendwhin eskaliert
                    else %Ziel erreicht
                        condi=0; %gehe aus der while schleife
                        %Abbruch532(j,LR)= Abbruch532(j,LR)+100;%disp('bringt nichts mehr'),
                    end
                    diffi = abs(tmp - BSR532mintrust);
                    if abs(diffi) < diffisoll, condi=0;  end  % normales Ende
                    
                    
                    
                    if hz532 >hzlim, condi = 0;
                        Abbruch532(j,LR)=Abbruch532(j,LR)+200; % es konvergiert nicht und landet immer bei zu hohen werten (1. Ansatz) --> Ermacht dann aber trotzdem weiter.
                    end   % Rdbeding. irrelevant
                    
                    if iter==itmax
                        Abbruch532(j,LR)=Abbruch532(j,LR)+100; %es konvergiert nicht (1. Ansatz) 
                    end
                    
                    % deltab
                    % qq=BSRAtFit532arr(j)
                end % 1. while
                
                
                if tmp < (BSR532mintrust-diffisoll)  % noetig falls iter > itmax. und er oben nicht auf ein ordentliches BSR gekommen ist
                    Btemp532=Btemp532./tmp.*BSR532mintrust;
                    %gesamtes Profil wird skaliert das UBC mindestens der dem
                    %mintrust entspricht
                end
                %wo532= find(Btemp532 < Wolkenschwelle532 );    %& H< hwo1);
                

                % alleaeroposi=sort(cat(1,wo532,wo532S,wo355));
                % da=diff(alleaeroposi);
                % da2=find(da >0);
                % woaerosol=[1, alleaeroposi(da2+1)];
                % stimmt das? Wir brauchen es nicht
                
                %Hier wird geguckt wo nur klare Luft ist und das LR f?r
                %klare Luft eingesetzt
                woaerosol = find(Btemp532 < Wolkenschwelle532);
                wowolke = find(Btemp532 > Wolkenschwelle532);
                LR532arr(woaerosol,j,LR) = LR532aerosol;
                %LR532Sarr(woaerosol,j) = LR532Saerosol;
                %LR355arr(woaerosol,j) = LR355aerosol;
                
                guteaeroposi= connrnge(H>  UeberlappEnde & H < hwo1 & Btemp532 < Wolkenschwelle532); % das sucht die zusammenhaengenden Bereiche ohne Wolken
                if length(guteaeroposi) > 1
                    guteaeroposi = (guteaeroposi(1):guteaeroposi(2));%geeingnete position fuer vergleich (kontrollrange fuer BSR355soll
                else
                    guteaeroposi = floor(UeberlappEnde/7.5 ):40;  % irgendwelche Positionen dicht unter Flugzeug
                    %ichmerkmirkomischepositionen = 1;
                    Abbruch532(j,LR) = Abbruch532(j,LR)+10; %keine Wolkenfreien Positionen gefunden, Einfach zwischen 200 und 300 m unterm Flugzeug zu CLEAR erkl?rt
                end
                
                
                % wir glauben, dass die Wolkenmaske und damit die Verteilung, an welchen
                % Positionen das LR f?r Aerosol oder Wolken gesetzt wurde, "richtig genug
                % ist" um die eigentliche Rechnung durchzuf?hren. Wir iterireren jetzt (mit
                % dem hoffentlich richtigem LR f?r Aerosol die Randbedingung so lange, bis
                % die L?sung f?r den Fall "keine dicke Wolke" stimmt.
                
                
                BSR532haben = mymedian(Btemp532(guteaeroposi));      %#ok<NASGU> % war mean  %% benutzt er das spaeter ueberhaupt noch?
                BSR532sollarr(j)=mymedian(BSR532Karlmedianvgl(guteaeroposi)); %BSRsoll ist jetzt das ziel fuer die upper boundary condition
                if BSR532sollarr(j) < BSR532mintrust-diffisoll  || ~isfinite(BSR532sollarr(j))
                    BSR532sollarr(j) = BSR532sollnotfall;
                    Abbruch532(j,LR) = Abbruch532(j,LR)+0.3; %????????????????????????????????????????????????das sollte gar nicht passieren k?nnen
                end
                
                
                
                condi=1; iter=0; itmax=500;
                controllBSRWert=zeros(itmax,1);
                while condi
                    iter=iter+1;
                    if iter >= itmax 
                        condi=0; 
                        disp('keine Konvergenz gefunden 532P');  j, 
                        Abbruch532(j,LR) = Abbruch532(j,LR)+10 %Konvergiert nicht 2.Ansatz
                        %hier w?re es jetzt gescchickt mit einem median
                        %CLidar aus vorherigen Messugnen Das signal zu
                        %multiplizieren und dann das als BSRsoll bzw. UBC
                        %zu nehmen
                    end
                    
                    BSRAtFiterr = BSRAtFit532arr(j) ./ 5;
                    [Beta, dBdR532, dBdLR532, dBdP532, CLidar532] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Klettnoise(Sel532P,j), ...
                        LR532arr(Sel532P,j,LR), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                    %Anmerkung, eingentlich ist C in der Lidargleichung nur
                    %ein fixer wert, aber weil wir die Lidargleichung f?r
                    %jedes rangegate umgestellt haben kommt ein Vektor
                    %raus. Da wo das komplett off ist stimmt mit dem
                    %retrieval etwas nicht. - z.B. H?ufig unterhalb des
                    %Bodens
                    
                    Betaaer532(Sel532P)=Beta-BeRa532(Sel532P,1);% nur der Aersosol anteil der Rueckstreuung, molekularer Anteil abgezogen
                    
                    Betaaer532err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Klettnoise(Sel532P,j));%Summe aller moeglichen Fehlerquellen
                    
                    %Fehler in BSR umgerechnet (warum auch immer)
                    Btemp532(Sel532P)=Beta./BeRa532(Sel532P,1); Btemp532err(Sel532P)=Betaaer532err(Sel532P,j)./BeRa532(Sel532P,1);
                    
                    BSR532haben = mymedian(Btemp532(guteaeroposi)); % war mymean % mittlerer Wert in UBC die dann mit mittlerem Wert mit veraenderter LBC verglichen wird
                    
                    q=BSRAtFit532arr(j)+0.2; %LBC veraendern
                    
                    [Beta2, ~, ~, ~, ~] = klettinv_ableit4( q, FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Klettnoise(Sel532P,j), ...
                        LR532arr(Sel532P,j,LR), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                    BetaAer532_2(Sel532P)=Beta2-BeRa532(Sel532P,1);
                    %fuer den Vergleich brauchen wir keine Fehlerabschaetzung
                    %BetaAer532_2err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Klettnoise(Sel532P,j));
                    Btemp532_2(Sel532P)=Beta2./BeRa532(Sel532P,1); %Btemp532_2err(Sel532P)=BetaAer532_2err(Sel532P,j)./BeRa532(Sel532P,1);
                    
                    BSR532haben2 = mymedian(Btemp532_2(guteaeroposi));   % war mean
                    %vergleich
                    deltab=(BSR532haben2-BSR532haben);
                    
                    %das ist damit man die zwischenwerte der iterationen
                    %plotten kann - eigentlich nicht wichtig
                    %                 BetaKlett2(:,iter,j)= Beta;
                    %                 LBCKlett2(:,j)=BSRAtFit532arr(j);
                    %                 UBCKlett2(:,iter,j)=  Btemp532;
                    %                 LRKlett2(:,iter,j)= LR532arr(Sel532P,j);
                    %                 CKlett2(:,iter,j)=CLidar532;
                    %                 clearKlett2(guteaeroposi,iter,j)= 1;
                    
                    HKlett2(:,j)=H(Sel532P); %#ok<AGROW,NASGU>
                    
                    
                    if abs(deltab) > 0.01 %wenn sich noch was veraendert hat
                        w=(BSR532sollarr(j) - BSR532haben) ./ deltab;
                        BSRAtFit532arr(j) = BSRAtFit532arr(j) + w.*0.2;
                    else
                        condi=0; %disp('bringt nichts mehr'),
                    end
                    
                    diffi = abs(BSR532haben - BSR532sollarr(j));
                    
                    if abs(diffi) < diffisoll, condi=0;  end  % normales Ende % das ist eigentlich doppelt gemoppelt, aber solange es funktioniert und nicht st?rt ... 
                    
                    %Das hier ist nur wenn man zwischendurch plotten will
                    %ob die ganzen iterationen ?berhaupt etwas bringen
                    controllBSRWert(iter) =  BSR532haben;
                end % while f?r die Randbedingung
                
                %sporadisch kommen zu niedrige Rd-bedingungen vor, weil vermutlich die
                %LR total falsch sind. Dies wird hier abgefangen: (gefunden bei 532S)
                %Wenn die untere Randbedingung winzig klein ist???
                if BSRAtFit532arr(j) < BSR532mintrust -diffisoll
                    BSRAtFit532arr(j) = BSR532mintrust;
                    Abbruch532(j,LR) = Abbruch532(j,LR)+2000;% ??? das passiert andauernd.  und vor allem macht er dann ja auch nichts mehr damit ??? vielleicht weil hier am Anfang sehr viel Wolke sehr weit oben ist? Oder ist die Wolke zu dicht und er interpretiert das Rauschen zu sehr.  
                end
                
                
                
                %LR variation hier rausgeloescht
                
                
                %matrizen werden zusammengesetzt
                %dritte Dimension 1 - 1.LR (theoretisches wasser), 2 -LR wie
                %Aerosol, 3 LR wie Eis in LAMPERT
                
                %molekulare anteil
                %BeRa532(Sel532P,1); % da das aus den monatsmitteln der Radiosonden kommt sollte das fuer jeden Flug gleich sein
                PositionUBC = zeros(size(BSR532Klett));
                
                %backscatterratio (beta total / beta mol)
                PositionUBC(guteaeroposi,j,LR)= 1;
                BSR532Klett(Sel532P,j,LR)=Btemp532(Sel532P);
                BSR532Kletterr(Sel532P,j,LR)=Btemp532err(Sel532P);
                
                %backscatter (mit oder ohne molekularem Anteil??? )
                BetaAer532Klett(Sel532P,j,LR)=Betaaer532(Sel532P);
                BetaAer532Kletterr(Sel532P,j,LR)=Betaaer532err(Sel532P,j);
                %attenuated backscatter
                attenu532P(Sel532P,j,LR) = P532Klett(Sel532P,j).*H(Sel532P).^2 ./ CLidar532(Sel532P(1));
                % Fehlerquelle im attenuated backscatter ist nur das Rauschen und ein fehler in der Annahme der oberen randbedingung
                attenu532Perr(Sel532P,j,LR) = abs(P532Klettnoise(Sel532P,j).*H(Sel532P).^2 ./ CLidar532(Sel532P(1))) + abs(attenu532P(Sel532P,j) ./ CLidar532(Sel532P(1)).*0.1.*CLidar532(Sel532P(1)));
                dBeta532dP(Sel532P,j,LR)=dBdP532; %Fehler durch Rauschen
                dBeta532dR(Sel532P,j,LR)=dBdR532; %Fehler durch untere Randbedingung
                dBeta532dLR(Sel532P,j,LR)=dBdLR532; %Fehler durch falsches LR
                %Fehler durch falsche obere Randbedingung ist linaer. wenn oben
                %10% mehr sind dann muesste das attenuated backscatter profil
                %komplett 10% mehr sein
                C532Lidar_mat(Sel532P,j,LR) =CLidar532; % LIdarkonstante, sollte theoretisch f?r alle drei LR gleich sein (Signal(rangekorrigiert und hintergrundbereinigt)*Lidarkonstante=attenuated backscatter. - guter sanity check
                Abbruch532_mat(j,LR)= Abbruch532(j);
                %attenuation
                AlphaAer532(:,j,LR)=BetaAer532Klett(:,j,LR).*LR532arr(:,j,LR);
                %attenuation ist Backscatter (beta)*LR
                AlphaAer532err(:,j,LR)= abs(BetaAer532Kletterr(:,j).*LR532arr(:,j,LR)) + abs(BetaAer532Klett(:,j).*LR532arrerr(:,j));
                
                
                
                
                
                % finales Abspeichern: ----------------------------------------
                
                %matrizen werden zusammengesetzt
                
               
                
                
                
                
                Wolkenmaske(wowolke,j,1) = 1;
                
                
            end % for Zeitschritte "entries"

        end % for LR
        
        
        
        P532roh=P532A'; %#ok<NASGU>
        P532Sroh=P532SA'; %#ok<NASGU>
        P355roh=P355A'; %#ok<NASGU>
        
        %das hier muss noch angepasst werden
        
        
        speichern532p=['BeRa532 AlRay532'...%rayleigh profile for molcular scattering
            ' P532roh P532Klett P532Klettnoise'... %rohsignal, Hintergrund korrigiertes signal, Noise
            ' BSRAtFit532arr LR532arr'... %Annahme untere Randbedingung (7.5-107.5m) Annahmen fuer Lidarratio
            ' C532Lidar_mat PositionUBC'...
            ' BetaAer532Klett BetaAer532Kletterr' ... %aerosol (and cloud) backscatter, basically non-molecular backscatter and its maximum error -namenvergessen-Fehlerfortpflanzung
            ' AlphaAer532 AlphaAer532err'...aerosol (and cloud) attenuation            % ' BSR532Klett BSR532Kletterr'... backscatter ratio (total backscatter / molcular backscatter  - unnoetig weil wir ja die einzelnen Komponenten davon speichern
            ' dBeta532dP  dBeta532dR dBeta532dLR'... %veraenderung von beta pro 1 veraenderung P (z.b. durch Rauschen), oder Randbedingung oder LR
            ' Abbruch532_mat P532Abackground']; %die Fehler variable, falls z.B. unterm Flugzeug eine Wolke ist oder es nicht konvergiert, die Hintergrundkorrektur, die wir abgekogen haben
        
        
        speichern532s=['P532Sroh P532SKlett P532SKlettnoise'];
        
        
        speichern355p=['BeRa355 AlRay355'...
            ' P355roh P355Klett P355Klettnoise'];
        
        sonstiges = [' H matlabzeit Wolkenmaske FitRangearr BSR532sollarr BSRAtFit532start' ];
        
        %         speichern532s=['P532SKlett P532SKletterr P532Sroh BSRAtFit532Sarr LR532Sarr BetaAer532SKlett BetaAer532SKletterr ' ...
        %             'AlphaAer532S AlphaAer532Serr BSR532SKlett BSR532SKletterr BSR532SKlettfest BSR532SKlettfesterr BetaAer532SKlettfest ' ...
        %             'BetaAer532SKlettfesterr dBeta532SdP  dBeta532SdR dBeta532SdLR Abbruch532S P532SAbackground BetaAer532SKlettausP BSR532SKlettausP'];
        %
        %
        %         speichern355p=['P355Klett P355Kletterr P355roh BSRAtFit355arr  LR355arr BetaAer355Klett BetaAer355Kletterr ' ...
        %             'AlphaAer355 AlphaAer355err BSR355Klett BSR355Kletterr BSR355Klettfest BSR355Klettfesterr BetaAer355Klettfest ' ...
        %             'BetaAer355Klettfesterr dBeta355dP  dBeta355dR dBeta355dLR Abbruch355 P355Abackground'];
        %
        
        
        
        liste=[speichern532p ' ' speichern532s ' ' speichern355p sonstiges];
        
        speichern =['save ' speicherdir Speichernamepraefix filename ' ' liste];
        
        eval(speichern)
        
        
    end % for alle Daten eines Tages
    
end % if files nicht leer



return



