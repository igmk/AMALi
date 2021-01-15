function ok = amali_eval_Wolkezumplotten(DatStr, aerofak, Speichernamepraefix, BSR532soll, BSRAtFit532start, Von, Bis)

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
% %%% campaign = 'ACLOUD' 'AFLUX' gibt es nicht mehr
% Speichernamepraefix = 'Tollertest'
% BSR532soll = 1.3;  angestrebter Referenzwert von BSR532 unter klaren
% Bedingungen   wenn BSR532soll = 'KARL' dann liest er den KARL-NyA Wert
% ein
% BSRAtFit532start   Anfangs(Rate)Wert der Randbedingung vom Klett
% Von / Bis: optional Die Datensatznummern in denen gerechnet werden soll
% wird Von / Bis nicht gesetzt rechnet er alles
% ok = amali_eval_Wolke('170523', 1, 'Test1', 1.3, 1.4, 2000, 2100)
% ok = amali_eval_Wolke('170523', 1, 'Test1', 'KARL', 1.4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ok=-1;

%
disp("Birtes Version")


% Einstellungen fuer Kampagnen
if strcmp(DatStr(1:2), '17')
    campaign = 'ACLOUD';
    Angstroem=1.4;
    speicherdir = '/atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2017/';
    amalidatendir='/lidar4/lidar/amali2011/data/mat/2017/' %#ok<NOPRT>
elseif strcmp(DatStr(1:2), '19')
    campaign = 'AFLUX';
    Angstroem=1.2;
    speicherdir = '/atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2019/';
    amalidatendir='/atm_meas/polar_5_6/amali/data/mat/2019/';
    
elseif strcmp(DatStr(1:2), '20')
    campaign = 'MOSAICACA';
    Angstroem=1.2;
    speicherdir = '/atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2020/';
    amalidatendir='/atm_meas/polar_5_6/amali/data/mat/2020/';
else
    disp('Campaign unknown')
    return
end
WvlExpo = 4-Angstroem;

if nargin > 5
    teilzeitjob = 1;
else
    teilzeitjob = 0;
end

% Definitionen ------------------------------------------------------------
Wvl532 =5.3207e-7;
Wvl355 =3.5471e-7;
LR532aerosol=23;
LR532wolke=15;
LR532Saerosol=23;
LR532Swolke=15;
LR355aerosol=23;
LR355wolke=15;
LRobergr = 120;        %
LRuntergr = 3;
LRerr = 2; % angenommener Fehler im LR  war 10
Wolkenschwelle532=5;
Wolkenschwelle532S=10;   % verrauschter, deshalb hoeher
Wolkenschwelle355=3;

tzlim = 5; hzlim=5;  % limit f?r Hoch- Tiefzahl bei Iterationenen, so oft versucht er es nochmal wenn LR auserhalb des erlaubten Berreiches liegt

% BSR at LBC
% BSRAtFit532start = 1.4;      % besser als Inputparameter
BSRAtFit355start = 1+(BSRAtFit532start-1) ./ 1.5.^WvlExpo; % rechnen 355er Randbedingung aus der von 532

%minimum BSR at UBC
BSR532mintrust  = 1.2;
BSR532Smintrust  = 1.2;
BSR355mintrust  = 1 + (BSR532mintrust -1) ./ 1.5.^WvlExpo;

% BSRAtFiterr=0.2; Fehler ist jetzt 20% von BSRAtFit
% BSR532soll = 1.3;  % ist Inputparameter
% BSR355soll kann erst berechnet werden, wenn BSR532(guteaeroposi) aus KARL
% bekannt ist
% BSR355soll = 1 + (BSR532soll-1) ./ 1.5.^WvlExpo;

% nur f?er den Notfall (NaN im KARL oder keine guteaeroposi) fuer UBC:
BSR532sollnotfall = 1.3;
BSR532Ssollnotfall = 1.3;
BSR355sollnotfall = 1 + (BSR532sollnotfall-1) ./ 1.5.^WvlExpo;


UeberlappEnde=300; % [m]
pretrigrange=1:400;
pretriggerbins=405;
Schwelle = 1e-8; % Signal (P) nur noch Rauschen

UeberlappEndeposi = round(UeberlappEnde ./ 7.5);

diffisoll = 0.05;  % so gut wollen wir BSR bestimmen
%FitRange=[2600, 2700]; %Standardwahl, wenn es nichts Besseres gibt

Hcalcrange =3500;    %muss hoch gesetzt werden wenn dsa Flugzeug ueber 3500 m fliegt



% ------------------------------------------------------------------------
% Rayleighstreuung aus median Radiosonde NyA, laedt meadian KARL
if strcmp(campaign, 'MOSAICACA')
    ptuinfile='/atm_meas/awipev/lidar/karl/matlab/ptu/2009.mat';
    load(ptuinfile)
    ozoinfile='/atm_meas/awipev/lidar/karl/matlab/ozo/2009.mat';
    load(ozoinfile)
    load /atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2020/aerosol_background_karl/KARLaverageBSR532ausSommer
    KarlH = H; % H wird sp?ter der AMALi Range-vektor
    BSR532Karlmean= BSR532mean; BSR532Karlmedian = BSR532median;
    BSR532SKarlmean= BSR532Smean; BSR532SKarlmedian = BSR532Smedian;
    BSR355Karlmean= BSR355mean; BSR355Karlmedian = BSR355median;
    load /atm_meas/polar_5_6/flight_data/gps/MOSAICACAGPS.txt -ascii
    flugzeit=MOSAICACAGPS(:,1);
    flughoehe=MOSAICACAGPS(:,2);
    wo=find(flughoehe > 10000);
    flughoehe(wo)=NaN;
elseif strcmp(campaign, 'ACLOUD')
    ptuinfile='/atm_meas/awipev/lidar/karl/matlab/ptu/1706.mat';
    load(ptuinfile)
    ozoinfile='/atm_meas/awipev/lidar/karl/matlab/ozo/1706.mat';
    load(ozoinfile)
    load /atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2017/aerosol_background_karl/KARLaverageBSR532ausACLOUD
    KarlH = H; % H wird sp?ter der AMALi Range-vektor
    BSR532Karlmean= BSR532mean; BSR532Karlmedian = BSR532median;
    BSR532SKarlmean= BSR532Smean; BSR532SKarlmedian = BSR532Smedian;
    BSR355Karlmean= BSR355mean; BSR355Karlmedian = BSR355median;
    load /atm_meas/polar_5_6/flight_data/gps/ACLOUDGPS.txt -ascii
    flugzeit=ACLOUDGPS(:,1);
    flughoehe=ACLOUDGPS(:,2);
    wo=find(flughoehe > 10000);
    flughoehe(wo)=NaN;
elseif strcmp(campaign, 'AFLUX')
    ptuinfile='/atm_meas/awipev/lidar/karl/matlab/ptu/1904.mat';
    load(ptuinfile)
    ozoinfile='/atm_meas/awipev/lidar/karl/matlab/ozo/1904.mat';
    load(ozoinfile)
    load /atm_meas/polar_5_6/amali/data/nadir_processed/cloud/2019/aerosol_background_karl/KARLaverageBSR532ausAFLUX
    KarlH = H; % H wird sp?ter der AMALi Range-vektor
    BSR532Karlmean= BSR532mean.*aerofak; BSR532Karlmedian = BSR532median.*aerofak;
    BSR532SKarlmean= BSR532Smean.*aerofak; BSR532SKarlmedian = BSR532Smedian.*aerofak;
    BSR355Karlmean= BSR355mean.*aerofak; BSR355Karlmedian = BSR355median.*aerofak;
    load /atm_meas/polar_5_6/flight_data/gps/AFLUXGPS.txt -ascii
    flugzeit=AFLUXGPS(:,1);
    flughoehe=AFLUXGPS(:,2);
    wo=find(flughoehe > 10000);
    flughoehe(wo)=NaN; %#ok<*FNDSB>
end

meanO3profile=(mymean(OZOO3Density'))';
Density = density(PTUPressure, PTUPressure./100, PTUTemperature,ones(size(PTUTemperature)).*2); %#ok<*NODEF>
%Extinktion nur Streuanteil (99% ausmacht)
%alpha ist gleich fuer beide polarisationen
PTUAlRay532 = Density .* raytotwq( Wvl532, PTUTemperature, Density);
PTUAlRay355 = Density .* raytotwq( Wvl355, PTUTemperature, Density);
% Rayleigh R?ckstreuung (die schwach polaris.- und temp- abhaengig ist
% wichtig: Streuung =  Luftdichte (Density) mal Wirkungsquerschnitt
PTUBeRa532=Density.*raybckwq (Wvl532,'p','p', PTUTemperature, Density);
PTUBeRa532S=Density.*raybckwq (Wvl532,'p','s', PTUTemperature, Density);
PTUBeRa355=Density.*raybckwq (Wvl355,'p','u', PTUTemperature, Density);

% Korrektur durch molek. Absorption, die praktisch nur durch O3
% vor allem wichig in Strato - hier vermutlich egal
ptudimen=size(PTUAlRay355);
for j=1:ptudimen(2)
    ARayab532=meanO3profile.*o3abswq(Wvl532,PTUTemperature(:,j));
    ARayab355=meanO3profile.*o3abswq(Wvl355,PTUTemperature(:,j));
    PTUAlRay532(:,j)=PTUAlRay532(:,j)+ARayab532;
    PTUAlRay355(:,j)=PTUAlRay355(:,j)+ARayab355;
end


% gemittelte Profile aus Ny-Alesund
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

tmp=real(log(NAlRay355));
wo=find(isnan(tmp(1:20)));
if ~isempty(wo)
    woend=wo(end);
    st=mymean(diff(tmp(woend+1:woend+11)));
    for jj=woend:-1:1
        tmp(jj)=tmp(woend+1)-(woend+1-jj).*st;
    end
end
NAlRay355 = exp(tmp);

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

tmp=real(log(NBeRa532S));
wo=find(isnan(tmp(1:20)));
if ~isempty(wo)
    woend=wo(end);
    st=mymean(diff(tmp(woend+1:woend+11)));
    for jj=woend:-1:1
        tmp(jj)=tmp(woend+1)-(woend+1-jj).*st;
    end
end
NBeRa532S = exp(tmp);


tmp=real(log(NBeRa355));
wo=find(isnan(tmp(1:20)));
if ~isempty(wo)
    woend=wo(end);
    st=mymean(diff(tmp(woend+1:woend+11)));
    for jj=woend:-1:1
        tmp(jj)=tmp(woend+1)-(woend+1-jj).*st;
    end
end
NBeRa355 = exp(tmp);



%-------------------------------------------------------------------------


% noch nicht fertig
% % Rayleighstreuung aus Dropsonde
% dropfilepath='/lidar4/bkulla/ACLOUD2017/dropsondes/2017/06/';
% day='/02/';
% dropname=''
% dropfilename=[dropfilepath day dropname]
% dropT=ncread(dropfilename,'');
% dropH=ncread(dropfilename,'');
% dropP=ncread(dropfilename,'');
% dropdens = density(dropP, dropP./100, dropT,ones(size(dropT)).*2);
%
% %Extinktion nur Streuanteil (der 99% ausmacht)
% AlRay532 = dropdens .* raytotwq( Wvl532, dropT, dropdens);
% AlRay355 = dropdens .* raytotwq( Wvl355, dropT, dropdens);
% % Rayleigh R?ckstreuung (die schwach polaris.- und temp- abhaengig ist
% % wichtig: Streuung =  Luftdichte (Density) mal Wirkungsquerschnitt
% BeRa532=dropdens.*raybckwq (Wvl532,'p','p', dropT, dropdens);
% BeRa532S=dropdens.*raybckwq (Wvl532,'p','s', dropT, dropdens);
% BeRa355=Density.*raybckwq (Wvl355,'p','p', PTUTemperature, Density);
% % Korrektur durch molek. Absorption, die praktisch nur durch O3
% % vor allem wichit in Strato - wenn es kein Ozon aus Dropsonde gibt
% % vernachl?ssigen wir diesen Quatsch
% %dropdimen=size(AlRay355);
% %for j=1:dropdimen(2)
% %ARayab532=meanO3profile.*o3abswq(Wvl532,PTUTemperature(:,j));
% %ARayab355=meanO3profile.*o3abswq(Wvl355,PTUTemperature(:,j));
% %AlRay532(:,j)=AlRay532(:,j)+ARayab532;
% %AlRay355(:,j)=AlRay355(:,j)+ARayab355;
% end
%



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
    % jetzt geht es aber mal so richtig los
    ok = NofFiles;
    for kk=1:NofFiles %Schleife ?ber alle Datens?tze
        disp(['laboriere an Datensatz ' num2str(kk) ' von ' num2str(ok)])
        load(files(kk,:))
        wostrich=findstr(files(kk,:),'/');
        if length(wostrich>1), wostrich=wostrich(end); end
        womat= findstr(files(kk,:), '.mat');
        filename = files(kk,wostrich+1:womat-1);
        datasize=size(alldata);
        entries=datasize(2);
        bins=datasize(3);
        data532=(reshape(alldata(1,:,:),entries, bins));   % Format umgedreht jetzt wie KARL
        data532s=(reshape(alldata(3,:,:),entries, bins));  % size =[H, Zeiten]
        data355=(reshape(alldata(5,:,:),entries, bins));
        data532c=(reshape(alldata(2,:,:),entries, bins));
        data532sc=(reshape(alldata(4,:,:),entries, bins));
        data355c=(reshape(alldata(6,:,:),entries, bins));
        
        
        matlabzeit = (allinfo(1,:)/86400) + datenum(1904,1,1);
        
        H=(0:7.5:3500)';    dH=7.5;    lim=50; %(Limit Signal hat Wolken)
        LH=length(H);
        P532A=zeros(entries,LH);
        P532SA=zeros(entries,LH);
        P355A=zeros(entries,LH);
        P532C=zeros(entries,LH);
        P532SC=zeros(entries,LH);
        P355C=zeros(entries,LH);
        
        % Im ersten Datensatz (k=1) wandeln wir die Karl H?hen und BSRWerte auf das
        % AMALi Gitter um
        if kk==1
            if BSR532soll=='KARL'
                [KarlHzumBoden,indexx]=sort(KarlH,'descend');
                BSR532KarlmedianzumBoden = BSR532Karlmedian(indexx);
                %BSR532KarlmeanzumBoden = BSR532Karlmean(indexx);
                BSR532SKarlmedianzumBoden = BSR532SKarlmedian(indexx);
                %BSR532SKarlmeanzumBoden = BSR532SKarlmean(indexx);
                BSR355KarlmedianzumBoden = BSR355Karlmedian(indexx);
                %BSR355KarlmeanzumBoden = BSR355Karlmean(indexx);
                % jetzt m?ssen wir das auf das H?hengitter des Amali interpolieren
                BSR532Karlmedianvgl = interp1(KarlHzumBoden, BSR532KarlmedianzumBoden, H);
                %BSR532Karlmeanvgl = interp1(KarlHzumBoden, BSR532KarlmeanzumBoden, H);
                BSR532SKarlmedianvgl = interp1(KarlHzumBoden, BSR532SKarlmedianzumBoden, H);
                %BSR532SKarlmeanvgl = interp1(KarlHzumBoden, BSR532SKarlmeanzumBoden, H);
                BSR355Karlmedianvgl = interp1(KarlHzumBoden, BSR355KarlmedianzumBoden, H);
                %BSR355Karlmeanvgl = interp1(KarlHzumBoden, BSR355KarlmeanzumBoden, H);
            else
                BSR532Karlmedianvgl =  BSR532sollnotfall;
                %BSR532Karlmeanvgl =  BSR532sollnotfall;
                BSR532SKarlmedianvgl =  BSR532Ssollnotfall;
                %BSR532SKarlmeanvgl =  BSR532Ssollnotfall;
                BSR355Karlmedianvgl =  BSR355sollnotfall;
                %BSR355Karlmeanvgl =  BSR355sollnotfall;
            end
        end
        
        P532Abackground=zeros(entries,1);
        P532SAbackground=zeros(entries,1);
        P355Abackground=zeros(entries,1);
        P532Cbackground=zeros(entries,1);
        P532SCbackground=zeros(entries,1);
        P355Cbackground=zeros(entries,1);
        
        
        %spatial offsets - personal communication Roland Neuber, Checked so
        %that cloud tops align. 
        vs532p=-1;
        vs532s=-1;
        vs355=+1;
        wo = find(data532c > 3.2e4);
        %data532creserve=data532c;
        data532c(wo)=NaN;
        wo = find(data532sc > 3.0e4);
        %data532screserve=data532sc;
        data532sc(wo)=NaN;
        wo = find(data355c > 3.2e4);
        %data355creserve=data355c;
        data355c(wo)=NaN;
        
        
        for j=1:entries
            P532Abackground(j) = mymean(data532(j,pretrigrange));
            P532A(j,:)=data532(j,pretriggerbins+1:pretriggerbins+LH)- P532Abackground(j);
            P532SAbackground(j) = mymean(data532s(j,pretrigrange));
            P532SA(j,:)=data532s(j,pretriggerbins+1:pretriggerbins+LH)-P532SAbackground(j);
            P355Abackground(j) = mymean(data355(j,pretrigrange));
            P355A(j,:)=data355(j,pretriggerbins+1:pretriggerbins+LH)- P355Abackground(j);
            P532C(j,:)=data532c(j,vs532p+pretriggerbins+1:vs532p+pretriggerbins+LH)-mymean(data532c(j,pretrigrange));
            P532SC(j,:)=data532sc(j,vs532s+pretriggerbins+1:vs532s+pretriggerbins+LH)-mymean(data532sc(j,pretrigrange));
            P355C(j,:)=data355c(j,vs355+pretriggerbins+1:vs355+pretriggerbins+LH)-mymean(data355c(j,pretrigrange));
            
        end
        
        %rechnen mit Lidarsigal PKlett --> transponiert, damit das format
        %ist wie bei Karl
        P532Klett =P532A';
        P532SKlett =P532SA';
        P355Klett= P355A';
        P532CKlett =P532C';
        P532SCKlett =P532SC';
        P355CKlett= P355C';
        P532fuerdepol= P532A';
        P532Sfuerdepol= P532SA';
        VolDep532 = P532Sfuerdepol ./ P532fuerdepol;
        
        %"gesamt" sind die zusammengesetzen Profile
        % diese in counts wie die pc-Profile
        wo532scompare=find(P532SCKlett < 10 & P532SCKlett > 4);
        fktoranapc532s = mymean(P532SKlett(wo532scompare))./mymean(P532SCKlett(wo532scompare));
        %P532Sklettgesamt = P532Sklett ./ fktoranapc532s;
        
        
        wo=find(P532Klett < Schwelle); P532Klett(wo)=Schwelle;
        wo=find(P532SKlett < Schwelle); P532SKlett(wo)=Schwelle;
        wo=find(P355Klett < Schwelle); P355Klett(wo)=Schwelle;
        
        P532Kletterr=zeros(size(P532Klett));
        backgroundnoise532=zeros(entries,1);
        P532SKletterr=zeros(size(P532SKlett));
        backgroundnoise532S=zeros(entries,1);
        P355Kletterr=zeros(size(P355Klett));
        backgroundnoise355=zeros(entries,1);
        
        
        for j=1:entries
            backgroundnoise532(j)=std(data532(j,pretrigrange)).*1.5;
            P532Kletterr(:,j)=real(sqrt(P532Klett(:,j)))+backgroundnoise532(j);
            backgroundnoise532S(j)=std(data532s(j,pretrigrange)).*1.5;
            P532SKletterr(:,j)=real(sqrt(P532SKlett(:,j)))+backgroundnoise532S(j);
            backgroundnoise355(j)=std(data355(j,pretrigrange)).*1.5;
            P355Kletterr(:,j)=real(sqrt(P355Klett(:,j)))+backgroundnoise355(j);
        end
        SNR532=P532Klett./P532Kletterr;
        SNR532S=P532SKlett./P532SKletterr;
        SNR355=P355Klett./P355Kletterr;
        
        
        % Rayleighgr?ssen + Density
        NAlRay532aufl=interp1(NH,NAlRay532,H);
        AlRay532=NAlRay532aufl(end:-1:1);
        NAlRay355aufl=interp1(NH,NAlRay355,H);
        AlRay355=NAlRay355aufl(end:-1:1);
        NBeRa532aufl=interp1(NH,NBeRa532,H);
        BeRa532=NBeRa532aufl(end:-1:1);
        NBeRa532Saufl=interp1(NH,NBeRa532S,H);
        BeRa532S=NBeRa532Saufl(end:-1:1);
        NBeRa355aufl=interp1(NH,NBeRa355,H);
        BeRa355=NBeRa355aufl(end:-1:1);
        NDensityaufl=interp1(NH,NDensity,H);
        Density=NDensityaufl(end:-1:1);
        
        
        % Nun der schlaue Klett
        % wir brauchen f?r jeden Kanal:
        % BSRAtFitarr = das array der BSR-Werte f?r den Fit
        % FitRange = die H?he von, bis in der im Durchschnitt BSR=BSRAtFitarr sein
        % soll z.B. [2000, 2400] % Upper boundary condition??? oder lower?? ist
        % dass hoehe oder distance?
        % H: Hoehenvektor
        % P: das Lidarprofil
        % Perr "error" Fehler des Lidarprofils
        % LRarr= array des Lidarverh?ltnisses, gleich gro? wie P
        % ALRay532 / BeRa532 die Rayleogh-Extinktions und R?ckstreuwerte
        % klettinv_ableit4 berechnet:
        % 1) Beta = Gesamtr?ckstreuung (aer+Rayleigh)
        % 2) dBdR die partielle Ableitung d Beta_aer / d Randbedingung (oben oder
        % unten???
        % 3) dBdLR= die part. Ableit. d Beta_aer / d Lidarverh?eltnis
        % 4) dBdP: die part. Abl. d Beta_aer / d Perr
        % 5) Vorschlag CLidar: die Lidarkonstane aus L?sung Beta, LR und elast.
        % Lidargleichung
        
        % Definitionen:
        % f?r alle Farben:
        FitRangearr=zeros(2,entries);
        Wolkenmaske=zeros(size(P532Klett)); %
        dimen=size(P532Klett);
        
        % f?r 532P
        LR532arr=ones(size(P532Klett)).*LR532wolke;    % es geht ja um Wolken
        LR532arrerr=ones(size(P532Klett)).*LRerr;
        BSRAtFit532arr=ones(entries,1).* BSRAtFit532start;
        Btemp532=zeros(LH,1); Btemp532_2=zeros(LH,1);    Btemp532err=Btemp532;
        Betaaer532=zeros(LH,1); Betaaer532err=Btemp532;  BetaAer532_2=zeros(LH,1);
        BetaAer532Klett=zeros(size(P532Klett)); BetaAer532Kletterr=zeros(size(P532Klett));
        BetaAer532Klettfest=zeros(size(P532Klett)); BetaAer532Klettfesterr=zeros(size(P532Klett));
        dBeta532dLR=zeros(size(P532Klett));
        dBeta532dR=zeros(size(P532Klett));
        dBeta532dP=zeros(size(P532Klett));
        BSR532Klett=zeros(size(P532Klett)); BSR532Kletterr=zeros(size(P532Klett));
        BSR532Klettfest=zeros(size(P532Klett)); BSR532Klettfesterr=zeros(size(P532Klett));
        C532Lidar=zeros(size(P532Klett));
        AlphaAer532= zeros(size(P532Klett));
        AlphaAer532err= zeros(size(P532Klett));
        attenu532P= zeros(size(P532Klett));       attenu532Perr= zeros(size(P532Klett));
        Abbruch532 = zeros(dimen(2),1);
        BSR532sollarr = zeros(dimen(2),1);
        % Abbruch532 gibt Abbruchkriterium an:
        % 0:hat gar nicht LR angepa?t
        % 1: normale Konvergenz
        % 2: LR = 5 (unteres Limit)
        % 3: LR = LRobergr (oberes Limit)
        % 4: hat sich totiteriert ohne Konvergenz
        % 5: guteaeroposi = 0;
        %
        bininlrschleife532=zeros(entries,1);
        Ptheo532=2.*BeRa532./H.^2.*exp(-2.*qdrupvar(H,2.*AlRay532)); %was passiert hier? Wo kommt die 2 her? ansonsten ist das die Lidargleichung
        
        
        % f?r 532S
        LR532Sarr=ones(size(P532Klett)).*LR532wolke;    % es geht ja um Wolken
        LR532Sarrerr=ones(size(P532Klett)).*LRerr;
        BSRAtFit532Sarr=ones(entries,1).* BSRAtFit532start;
        Btemp532S=zeros(LH,1); Btemp532S_2=zeros(LH,1);    Btemp532Serr=Btemp532S;
        Betaaer532S=zeros(LH,1); Betaaer532Serr=Btemp532S;  BetaAer532S_2=zeros(LH,1);
        BetaAer532SKlett=zeros(size(P532SKlett)); BetaAer532SKletterr=zeros(size(P532SKlett));
        BetaAer532SKlettfest=zeros(size(P532SKlett)); BetaAer532SKlettfesterr=zeros(size(P532SKlett));
        BetaAer532SKlettausP=zeros(size(P532SKlett)); BetaAer532SKlettausPerr=zeros(size(P532SKlett));
        dBeta532SdLR=zeros(size(P532SKlett));
        dBeta532SdR=zeros(size(P532SKlett));
        dBeta532SdP=zeros(size(P532SKlett));
        BSR532SKlett=zeros(size(P532SKlett)); BSR532SKletterr=zeros(size(P532SKlett));
        BSR532SKlettfest=zeros(size(P532SKlett)); BSR532SKlettfesterr=zeros(size(P532SKlett));
        BSR532SKlettausP=zeros(size(P532SKlett)); BSR532SKlettausPerr=zeros(size(P532SKlett));
        C532SLidar=zeros(size(P532SKlett));
        AlphaAer532S= zeros(size(P532SKlett));
        AlphaAer532Serr= zeros(size(P532SKlett));
        attenu532S= zeros(size(P532SKlett));       attenu532Serr= zeros(size(P532SKlett));
        Abbruch532S = zeros(dimen(2),1);
        BSR532Ssollarr = zeros(dimen(2),1);
        % Abbruch532 gibt Abbruchkriterium an:
        % 0:hat gar nicht LR angepa?t
        % 1: normale Konvergenz
        % 2: LR = 5 (unteres Limit)
        % 3: LR = LRobergr (oberes Limit)
        % 4: hat sich totiteriert ohne Konvergenz
        % 5: guteaeroposi = 0;
        %
        bininlrschleife532S=zeros(entries,1);
        Ptheo532S=2.*BeRa532S./H.^2.*exp(-2.*qdrupvar(H,2.*AlRay532));
        
        
        % f?r 355P
        LR355arr=ones(size(P355Klett)).*LR355wolke;    % es geht ja um Wolken
        LR355arrerr=ones(size(P355Klett)).*LRerr;
        BSRAtFit355arr=ones(entries,1).* BSRAtFit355start;
        Btemp355=zeros(LH,1); Btemp355_2=zeros(LH,1);    Btemp355err=Btemp355;
        Betaaer355=zeros(LH,1); Betaaer355err=Btemp355;  BetaAer355_2=zeros(LH,1);
        BetaAer355Klett=zeros(size(P355Klett)); BetaAer355Kletterr=zeros(size(P355Klett));
        BetaAer355Klettfest=zeros(size(P355Klett)); BetaAer355Klettfesterr=zeros(size(P355Klett));
        dBeta355dLR=zeros(size(P355Klett));
        dBeta355dR=zeros(size(P355Klett));
        dBeta355dP=zeros(size(P355Klett));
        BSR355Klett=zeros(size(P355Klett)); BSR355Kletterr=zeros(size(P355Klett));
        BSR355Klettfest=zeros(size(P355Klett)); BSR355Klettfesterr=zeros(size(P355Klett));
        C355Lidar=zeros(size(P355Klett));
        AlphaAer355= zeros(size(P355Klett));
        AlphaAer355err= zeros(size(P355Klett));
        attenu355P= zeros(size(P355Klett));       attenu355Perr= zeros(size(P355Klett));
        Abbruch355 = zeros(dimen(2),1);
        BSR355sollarr = zeros(dimen(2),1);
        % Abbruch355 gibt Abbruchkriterium an:
        % 0:hat gar nicht LR angepa?t
        % 1: normale Konvergenz
        % 2: LR = 5 (unteres Limit)
        % 3: LR = LRobergr (oberes Limit)
        % 4: hat sich totiteriert ohne Konvergenz
        % 5: guteaeroposi = 0;
        %
        bininlrschleife355=zeros(entries,1);
        Ptheo355=2.*BeRa355./H.^2.*exp(-2.*qdrupvar(H,2.*AlRay355));
        
        
        
        if teilzeitjob==1
            start = Von; ende = Bis;
        else
            start = 1; ende = entries;
        end
        teiler=floor(ende./20);
        
        
        
        % Hier geht es los: die gro?e Schleife ueber alle Zeitschritte
        %
        %
        
        for j=start:ende   %1:entries    %2000:5000  start:ende
            
            if (mod(j,teiler)) ==0, disp(['Profile: ' num2str(j) ' von ' num2str(ende)]), end
            if j ==5000
                q=1;
            end
            
            
            % Definitionen innerhalb eines Zeitschrittes
            Psoll532=Ptheo532 ./Ptheo532(UeberlappEndeposi).*P532Klett(UeberlappEndeposi,j);
            Psoll532S=Ptheo532S ./Ptheo532S(UeberlappEndeposi).*P532SKlett(UeberlappEndeposi,j);
            Psoll355=Ptheo355 ./Ptheo355(UeberlappEndeposi).*P355Klett(UeberlappEndeposi,j);
            ichmerkmirkomischepositionen = 0;
            
            
            
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
            
            % "select" die guten H?henbins, in denen gerechnet werden soll.
            %
            Sel532S = connrnge( H >= 0 & H <= Hcalcrange & ...
                P532SKlett(:,j) > 0 & ...
                Density(:,1) > 0, 1);
            if length(Sel532S) > 1
                Sel532S = (Sel532S(1):Sel532S(2));
            else
                Sel532S = [];
            end
            
            
            % "select" die guten H?henbins, in denen gerechnet werden soll.
            %
            Sel355P = connrnge( H >= 0 & H <= Hcalcrange & ...
                P355Klett(:,j) > 0 & ...
                Density(:,1) > 0, 1);
            if length(Sel355P) > 1
                Sel355P = (Sel355P(1):Sel355P(2));
            else
                Sel355P = [];
            end
            
            
            
            
            
            
            
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
            [VoI,BiI]= minintsuche(abs(Psoll532(Sel532P(40:200))-P532Klett(Sel532P(40:200),j)),40);
            clearint=VoI:BiI;
            while condi
                iter = iter+1;
                if iter > itmax, condi = 0; end
                
                BSRAtFiterr = BSRAtFit532arr(j) ./ 5; %warum wird hier durch 5 geteilt irgendwas fuer fehlerabschaetzung
                [Beta, dBdR532, dBdLR532, dBdP532, CLidar] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Kletterr(Sel532P,j), ...
                    LR532arr(Sel532P,j), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                Betaaer532(Sel532P)=Beta-BeRa532(Sel532P,1);
                Betaaer532err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Kletterr(Sel532P,j));
                Btemp532(Sel532P)=Beta./BeRa532(Sel532P,1); Btemp532err(Sel532P)=Betaaer532err(Sel532P,j)./BeRa532(Sel532P,1);
                tmp=mymean(Btemp532(Sel532P(clearint)));
                
                [Beta2, ~, ~, ~, ~] = klettinv_ableit4( BSRAtFit532arr(j)+1, FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Kletterr(Sel532P,j), ...
                    LR532arr(Sel532P,j), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                BetaAer532_2(Sel532P)=Beta2-BeRa532(Sel532P,1);
                %BetaAer532_2err(Sel532P,j)=abs(dBdR5322.*BSRAtFiterr)+abs(dBdLR5322.*LR532arrerr(Sel532P,j))+abs(dBdP5322.*P532Kletterr(Sel532P,j));
                Btemp532_2(Sel532P)=Beta2./BeRa532(Sel532P,1); %Btemp532_2err(Sel532P)=BetaAer532_2err(Sel532P,j)./BeRa532(Sel532P,1);
                tmp2=mymean(Btemp532_2(Sel532P(clearint)));
                
                deltasoll = BSR532mintrust - tmp;
                deltab=(tmp2-tmp);
                
                %Zum plotten der zwischenschritte
                 BetaKlett1(:,iter,j)= Beta;
                 LBCKlett1(:,j)=BSRAtFit532arr(j);
                 UBCKlett1(:,iter,j)=  Btemp532;
                 CKlett1(:,iter,j)=CLidar;
                 clearKlett1(:,iter,j)=clearint;
                 HKlett1(:,j)=H(Sel532P);
                 
                if abs(deltasoll) > 0.01
                    %w=(BSR532sollarr(j) - tmp) ./ deltab;
                    BSRAtFit532arr(j) = BSRAtFit532arr(j) +1./deltab.*deltasoll;
                    if BSRAtFit532arr(j) > 1e6, hz532=hz532+1; BSRAtFit532arr(j) =1e6; end
                else
                    condi=0; 
                    Abbruch532(j)= Abbruch532(j)+100;%disp('bringt nichts mehr'),
                end
                diffi = abs(tmp - BSR532mintrust);
                if abs(diffi) < diffisoll, condi=0;  end  % normales Ende
                if hz532 >hzlim, condi = 0;
                    Abbruch532(j)=Abbruch532(j)+200
                end   % Rdbeding. irrelevant
                
                % deltab
                % qq=BSRAtFit532arr(j)
            end % 1. while
            
            
            if tmp < (BSR532mintrust-diffisoll)  % noetig falls iter > itmax.
                Btemp532=Btemp532./tmp.*BSR532mintrust;
                %gesamtes Profil wird skaliert das UBC mindestens der dem
                %mintrust entspricht
            end
            wo532= find(Btemp532 < Wolkenschwelle532 );    %& H< hwo1);
      
            
            
            % alleaeroposi=sort(cat(1,wo532,wo532S,wo355));
            % da=diff(alleaeroposi);
            % da2=find(da >0);
            % woaerosol=[1, alleaeroposi(da2+1)];
            % stimmt das? Wir brauchen es nicht
            
            woaerosol = find(Btemp532 < Wolkenschwelle532);
            LR532arr(woaerosol,j) = LR532aerosol;
            LR532Sarr(woaerosol,j) = LR532Saerosol;
            LR355arr(woaerosol,j) = LR355aerosol;
            guteaeroposi= connrnge(H>  UeberlappEnde & H < hwo1 & Btemp532 < Wolkenschwelle532); % das sucht den laengsten zusammenhaengenden Bereich ohne Wolken
            if length(guteaeroposi) > 1
                guteaeroposi = (guteaeroposi(1):guteaeroposi(2));%geeingnete position fuer vergleich (kontrollrange fuer BSR355soll
            else
                guteaeroposi = 20:40;  % irgendwelche Positionen dicht unter Flugzeug
                ichmerkmirkomischepositionen = 1;
                Abbruch532 = Abbruch532+10
            end
            
            
            % wir glauben, dass die Wolkenmaske und damit die Verteilung, an welchen
            % Positionen das LR f?r Aerosol oder Wolken gesetzt wurde, "richtig genug
            % ist" um die eigentliche Rechnung durchzuf?hren. Wir iterireren jetzt (mit
            % dem hoffentlich richtigem LR f?r Aerosol die Randbedingung so lange, bis
            % die L?sung f?r den Fall "keine dicke Wolke" stimmt.
            
            
            BSR532haben = mymedian(Btemp532(guteaeroposi));      % war mean
            BSR532sollarr(j)=mymedian(BSR532Karlmedianvgl(guteaeroposi)); %BSRsoll ist jetzt das ziel fuer die upper boundary condition
            if BSR532sollarr(j) < BSR532mintrust-diffisoll  || ~isfinite(BSR532sollarr(j))
                BSR532sollarr(j) = BSR532sollnotfall;
                Abbruch532(j) = Abbruch532(j)+0.3; %????????????????????????????????????????????????
            end
            
            
            condi=1; iter=0; itmax=500;
            controllBSRWert=zeros(itmax,1);
            while condi
                iter=iter+1;
                if iter >= itmax, condi=0; disp('trotz langer, muehsamer Suche keine Konvergenz gefunden 532P');  j, end
                
                BSRAtFiterr = BSRAtFit532arr(j) ./ 5;
                [Beta, dBdR532, dBdLR532, dBdP532, CLidar] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Kletterr(Sel532P,j), ...
                    LR532arr(Sel532P,j), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                Betaaer532(Sel532P)=Beta-BeRa532(Sel532P,1);
                Betaaer532err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Kletterr(Sel532P,j));
                Btemp532(Sel532P)=Beta./BeRa532(Sel532P,1); Btemp532err(Sel532P)=Betaaer532err(Sel532P,j)./BeRa532(Sel532P,1);
                BSR532haben = mymedian(Btemp532(guteaeroposi)); % war mymean
                
                q=BSRAtFit532arr(j)+0.2;
                
                [Beta2, dBdR532, dBdLR532, dBdP532, CLidar] = klettinv_ableit4( q, FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Kletterr(Sel532P,j), ...
                    LR532arr(Sel532P,j), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                BetaAer532_2(Sel532P)=Beta2-BeRa532(Sel532P,1);
                BetaAer532_2err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Kletterr(Sel532P,j));
                Btemp532_2(Sel532P)=Beta2./BeRa532(Sel532P,1); Btemp532_2err(Sel532P)=BetaAer532_2err(Sel532P,j)./BeRa532(Sel532P,1);
                BSR532haben2 = mymedian(Btemp532_2(guteaeroposi));   % war mean
                deltab=(BSR532haben2-BSR532haben);
                
                
                 BetaKlett2(:,iter,j)= Beta;
                 LBCKlett2(:,j)=BSRAtFit532arr(j);
                 UBCKlett2(:,iter,j)=  Btemp532;
                 LRKlett2(:,iter,j)= LR532arr(Sel532P,j);
                 CKlett2(:,iter,j)=CLidar;
                 clearKlett2(:,iter,j)=guteaeroposi;
                 HKlett2(:,j)=H(Sel532P);
                 
                 
                if abs(deltab) > 0.01
                    w=(BSR532sollarr(j) - BSR532haben) ./ deltab;
                    BSRAtFit532arr(j) = BSRAtFit532arr(j) + w.*0.2;
                else
                    condi=0; %disp('bringt nichts mehr'),
                end
                diffi = abs(BSR532haben - BSR532sollarr(j));
                if abs(diffi) < diffisoll, condi=0;  end  % normales Ende
                
                controllBSRWert(iter) =  BSR532haben;
            end % while f?r die Randbedingung
            
            %sporadisch kommen zu niedrige Rd-bedingungen vor, weil vermutlich die
            %LR total falsch sind. Dies wird hier abgefangen: (gefunden bei 532S)
            if BSRAtFit532arr(j) < BSR532mintrust -diffisoll
                BSRAtFit532arr(j) = BSR532mintrust;
                Abbruch532(j) = Abbruch532(j)+20;
            end
            
            
            
            
            % Das ist jetzt der Casus Knacktus, Iteration des LR
            % wahrscheinlich ist beta immer noch zu klein, denn bei dicken Wolken
            % hat der Wert BSRATFit keinen Einflu? auf die L?sung
            % Das bedeutet, die einzige Chance, eine realistische L?sung zu
            % bekommen besteht darin, iterativ das LR in der Wolke solange
            % anzupassen, bis die L?sung stimmt. Dies geschieht hier
            
            %das ist sinnlos. wir sagen oben soll es stimmen, aber oben ist
            %vom LR ueberhaupt nicht beeinflusst. 
            
            
            
            % wir passen vor der finalen Iteration die Wolkenmaske noch einmal an
            wowolke = find(Btemp532 > Wolkenschwelle532);
            if ~isempty(wowolke)% sonst ja sinnlos
                
                woaerosol = find(Btemp532 <= Wolkenschwelle532 & Btemp532S <= Wolkenschwelle532S & Btemp355 <= Wolkenschwelle355);
                LR532arr(wowolke,j)= LR532wolke;
                LR532arr(woaerosol,j)= LR532aerosol;
                LR532Sarr(wowolke,j)= LR532Swolke;
                LR532Sarr(woaerosol,j)= LR532Saerosol;
                LR355arr(wowolke,j)= LR355wolke;
                LR355arr(woaerosol,j)= LR355aerosol;
                
                % 532P
                condi=1; iter=0; itmax=500;
                hz532=0; tz532=0;  % hz532 / tz532 Zaehler fuer hohe - tiefe Werte
                gut = 1; % Prozent-Abweichung wie gut muss BSR an BSR532soll kommen
                refer=zeros(itmax,1);  refer2=zeros(itmax,1);
                while condi
                    iter = iter+1;
                    if iter ==1, bininlrschleife532(j) = 1; end
                    if iter > 500, condi=0; disp('keine Konvergenz!'); Abbruch532(j)=Abbruch532(j)+4;   end
                    
                    [Beta, dBdR532, dBdLR532, dBdP532, CLidar532] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Kletterr(Sel532P,j), ...
                        LR532arr(Sel532P,j), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                    Betaaer532(Sel532P)=Beta-BeRa532(Sel532P,1);
                    Betaaer532err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Kletterr(Sel532P,j));
                    Btemp532(Sel532P)=Beta./BeRa532(Sel532P,1); Btemp532err(Sel532P)=Betaaer532err(Sel532P,j)./BeRa532(Sel532P,1);
                    BSR532haben = mymedian(Btemp532(guteaeroposi));      %% war mean
                    
                    LRtest=LR532arr(:,j); LRtest(wowolke)=LRtest(wowolke)+1;
                    
                    [Beta2, dB2dR, dB2dLR, dB2dP, CLidar2] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Kletterr(Sel532P,j), ...
                        LRtest(Sel532P), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
                    BetaAer532_2(Sel532P)=Beta2-BeRa532(Sel532P,1);
                    BetaAer532_2err(Sel532P,j)=abs(dB2dR.*BSRAtFiterr)+abs(dB2dLR.*LR532arrerr(Sel532P,j))+abs(dB2dP.*P532Kletterr(Sel532P,j));
                    Btemp532_2(Sel532P)=Beta2./BeRa532(Sel532P,1); Btemp532_2err(Sel532P)=BetaAer532_2err(Sel532P,j)./BeRa532(Sel532P,1);
                    BSR532haben2 = mymedian(Btemp532_2(guteaeroposi));    % mean
                    
                    refer(iter) = BSR532haben;
                    refer2(iter) = BSR532haben2;
                    dxdlr = refer2(iter)-refer(iter);
                    diffi = BSR532sollarr(j) - refer(iter);
                    deltax=diffi./dxdlr;
                    
                 BetaKlett3(:,iter,j)= Beta;
                 LBCKlett3(:,j)=BSRAtFit532arr(j);
                 UBCKlett3(:,iter,j)=  Btemp532;
                 LRKlett3(:,iter,j)= LR532arr(Sel532P,j);
                 CKlett3(:,iter,j)=CLidar;
                 clearKlett3(:,iter,j)=guteaeroposi;
                 HKlett3(:,j)=H(Sel532P);
                 
                 
                    if refer(iter) ./ BSR532sollarr(j) > 1+gut ./ 100          % BSR zu gro?
                        LR532arr(wowolke,j) = LR532arr(wowolke,j) +deltax;
                        if mymean(LR532arr(wowolke,j)) > LRobergr
                            LR532arr(wowolke,j) = LRobergr; hz532=hz532+1;
                            if (BSRAtFit532arr(j) > 2 && hz532 < hzlim), BSRAtFit532arr(j) = BSRAtFit532arr(j)./1.1; end
                        end
                        if mymean(LR532arr(wowolke,j)) < LRuntergr
                            LR532arr(wowolke,j) = 5; tz532=tz532+1;
                            if (BSRAtFit532arr(j) < 1e6 && tz532 < tzlim), BSRAtFit532arr(j) = BSRAtFit532arr(j).*1.1; end
                        end
                    elseif refer(iter) ./ BSR532sollarr(j) < 1-gut ./ 100      % BSR zu klein
                        LR532arr(wowolke,j) = LR532arr(wowolke,j) +deltax;
                        if mymean(LR532arr(wowolke,j)) > LRobergr
                            LR532arr(wowolke,j) = LRobergr; hz532=hz532+1;
                            if (BSRAtFit532arr(j) > 2 && hz532 < hzlim), BSRAtFit532arr(j) = BSRAtFit532arr(j)./1.1; end
                        end
                        if mymean(LR532arr(wowolke,j)) < LRuntergr
                            LR532arr(wowolke,j) = 5; tz532=tz532+1;
                            if (BSRAtFit532arr(j) < 1e6 && tz532 < tzlim), BSRAtFit532arr(j) = BSRAtFit532arr(j).*1.1; end
                        end
                    else
                        condi=0; Abbruch532(j)=Abbruch532(j)+1;  % das sch?ne Ende
                    end
                    
                    
                    if hz532>=hzlim, condi = 0; Abbruch532(j)=Abbruch532(j)+3;  end      % disp('hz > 3'); i, end
                    if tz532>=tzlim, condi = 0; Abbruch532(j)=Abbruch532(j)+2;  end      % disp('tz > 3'); i, end
                    
                end    % while des LR
                
                
               
                
                
                
            end % if wowolke nicht leer
            
            
            if (ichmerkmirkomischepositionen ~=0)
                Abbruch532(j)=Abbruch532(j)+5;
                Abbruch532S(j)=Abbruch532S(j)+5;
                Abbruch355(j)=Abbruch355(j)+5;
            end
            
            % Wir reduzieren den Fehler des LR (LRerr) wenn Abbruch = 1 ist
            if Abbruch532(j) ==1
                LR532arrerr(wowolke,j) = 1;
                Betaaer532err(Sel532P,j)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Kletterr(Sel532P,j));
            end
            
            if Abbruch532S(j) ==1
                LR532Sarrerr(wowolke,j)=1;
                Betaaer532Serr(Sel532S,j)=abs(dBdR532S.*BSRAtFiterr)+abs(dBdLR532S.*LR532Sarrerr(Sel532S,j))+abs(dBdP532S.*P532SKletterr(Sel532S,j));
            end
            
            
            if Abbruch355(j) ==1
                LR355arrerr(wowolke,j)=1;
                Betaaer355err(Sel355P,j)=abs(dBdR355.*BSRAtFiterr)+abs(dBdLR355.*LR355arrerr(Sel355P,j))+abs(dBdP355.*P355Kletterr(Sel355P,j));
            end
            
            
            
            % finales Abspeichern:
            
            BSR532Klett(Sel532P,j)=Btemp532(Sel532P); BSR532Kletterr(Sel532P,j)=Btemp532err(Sel532P);
            BetaAer532Klett(Sel532P,j)=Betaaer532(Sel532P); BetaAer532Kletterr(Sel532P,j)=Betaaer532err(Sel532P,j);
            attenu532P(Sel532P,j) = P532Klett(Sel532P,j).*H(Sel532P).^2 ./ CLidar532(Sel532P(1));
            attenu532Perr(Sel532P,j) = abs(P532Kletterr(Sel532P,j).*H(Sel532P).^2 ./ CLidar532(Sel532P(1))) + abs(attenu532P(Sel532P,j) ./ CLidar532(Sel532P(1)).*0.1.*CLidar532(Sel532P(1)));
            dBeta532dP(Sel532P,j)=dBdP532;
            dBeta532dR(Sel532P,j)=dBdR532;
            dBeta532dLR(Sel532P,j)=dBdLR532;
            C532Lidar(Sel532P,j) =CLidar532;
            AlphaAer532(:,j)=BetaAer532Klett(:,j).*LR532arr(:,j);
            AlphaAer532err(:,j)= abs(BetaAer532Kletterr(:,j).*LR532arr(:,j)) + abs(BetaAer532Klett(:,j).*LR532arrerr(:,j));
            
            BSR532SKlett(Sel532S,j)=Btemp532S(Sel532S); BSR532SKletterr(Sel532S,j)=Btemp532Serr(Sel532S);
            BetaAer532SKlett(Sel532S,j)=Betaaer532S(Sel532S); BetaAer532SKletterr(Sel532S,j)=Betaaer532Serr(Sel532S,j);
            attenu532S(Sel532S,j) = P532SKlett(Sel532S,j).*H(Sel532S).^2 ./ CLidar532S(Sel532S(1));
            attenu532Serr(Sel532S,j) = abs(P532SKletterr(Sel532S,j).*H(Sel532S).^2 ./ CLidar532S(Sel532S(1))) + abs(attenu532S(Sel532S,j) ./ CLidar532S(Sel532S(1)).*0.1.*CLidar532S(Sel532S(1)));
            dBeta532SdP(Sel532S,j)=dBdP532S;
            dBeta532SdR(Sel532S,j)=dBdR532S;
            dBeta532SdLR(Sel532S,j)=dBdLR532S;
            C532SLidar(Sel532S,j) =CLidar532S;
            AlphaAer532S(:,j)=BetaAer532SKlett(:,j).*LR532Sarr(:,j);
            AlphaAer532Serr(:,j)= abs(BetaAer532SKletterr(:,j).*LR532Sarr(:,j)) + abs(BetaAer532SKlett(:,j).*LR532Sarrerr(:,j));
            
            BSR355Klett(Sel355P,j)=Btemp355(Sel355P); BSR355Kletterr(Sel355P,j)=Btemp355err(Sel355P);
            BetaAer355Klett(Sel355P,j)=Betaaer355(Sel355P); BetaAer355Kletterr(Sel355P,j)=Betaaer355err(Sel355P,j);
            attenu355P(Sel355P,j) = P355Klett(Sel355P,j).*H(Sel355P).^2 ./ CLidar355(Sel355P(1));
            attenu355Perr(Sel355P,j) = abs(P355Kletterr(Sel355P,j).*H(Sel355P).^2 ./ CLidar355(Sel355P(1))) + abs(attenu355P(Sel355P,j) ./ CLidar355(Sel355P(1)).*0.1.*CLidar355(Sel355P(1)));
            dBeta355dP(Sel355P,j)=dBdP355;
            dBeta355dR(Sel355P,j)=dBdR355;
            dBeta355dLR(Sel355P,j)=dBdLR355;
            C355Lidar(Sel355P,j) =CLidar355;
            AlphaAer355(:,j)=BetaAer355Klett(:,j).*LR355arr(:,j);
            AlphaAer355err(:,j)= abs(BetaAer355Kletterr(:,j).*LR355arr(:,j)) + abs(BetaAer355Klett(:,j).*LR355arrerr(:,j));
            
            
            
            
            
            Wolkenmaske(wowolke,j) = 1;
            
            % egal ob Wolke oder nicht - jetzt rechnen wir nochmal mit LR532wolke und
            % der Rd-beding, die ggf. gefunden wurde f?r die "fest" L?sungen
            %Warum???? --- vor allem ueber schreibt das alles.
            
            LRfest = ones(size(LR532arr)).*LR532wolke;
            [Beta, dBdR532, dBdLR532, dBdP532, CLidar] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532P), P532Klett(Sel532P,j), P532Kletterr(Sel532P,j), ...
                LRfest(Sel532P,j), AlRay532(Sel532P,1), BeRa532(Sel532P,1));
            Betaaer532(Sel532P)=Beta-BeRa532(Sel532P,1);
            Betaaer532err(Sel532P,1)=abs(dBdR532.*BSRAtFiterr)+abs(dBdLR532.*LR532arrerr(Sel532P,j))+abs(dBdP532.*P532Kletterr(Sel532P,j));
            Btemp532(Sel532P)=Beta./BeRa532(Sel532P,1); Btemp532err(Sel532P)=Betaaer532err(Sel532P,1)./BeRa532(Sel532P,1);
            BSR532Klettfest(Sel532P,j)=Btemp532(Sel532P); BSR532Klettfesterr(Sel532P,j)=Btemp532err(Sel532P);
            BetaAer532Klettfest(Sel532P,j)=Betaaer532(Sel532P); BetaAer532Klettfesterr(Sel532P,j)=Betaaer532err(Sel532P,1);
            
            LRfest = ones(size(LR532Sarr)).*LR532wolke;
            [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532arr(j), FitRangearr(:,j), H(Sel532S), P532SKlett(Sel532S,j), P532SKletterr(Sel532S,j), ...
                LRfest(Sel532S,j), AlRay532(Sel532S,1), BeRa532S(Sel532S,1));
            Betaaer(Sel532S)=Beta-BeRa532S(Sel532S,1);
            Betaaererr(Sel532S,1)=abs(dBdR.*BSRAtFiterr)+abs(dBdLR.*LR532Sarrerr(Sel532S,j))+abs(dBdP.*P532Kletterr(Sel532S,j));
            Btemp(Sel532S)=Beta./BeRa532S(Sel532S,1); Btemperr(Sel532S)=Betaaererr(Sel532S,1)./BeRa532S(Sel532S,1);
            BSR532SKlettfest(Sel532S,j)=Btemp(Sel532S); BSR532SKlettfesterr(Sel532S,j)=Btemperr(Sel532S);
            BetaAer532SKlettfest(Sel532S,j)=Betaaer(Sel532S); BetaAer532SKlettfesterr(Sel532S,j)=Betaaererr(Sel532S,1);
            
            
            LRfest = ones(size(LR355arr)).*LR355wolke;
            [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit355arr(j), FitRangearr(:,j), H(Sel355P), P355Klett(Sel355P,j), P355Kletterr(Sel355P,j), ...
                LRfest(Sel355P,j), AlRay355(Sel355P,1), BeRa355(Sel355P,1));
            Betaaer(Sel355P)=Beta-BeRa355(Sel355P,1);
            Betaaererr(Sel355P,1)=abs(dBdR.*BSRAtFiterr)+abs(dBdLR.*LR355arrerr(Sel355P,j))+abs(dBdP.*P355Kletterr(Sel355P,j));
            Btemp(Sel355P)=Beta./BeRa355(Sel355P,1); Btemperr(Sel355P)=Betaaererr(Sel355P,1)./BeRa355(Sel355P,1);
            BSR355Klettfest(Sel355P,j)=Btemp(Sel355P); BSR355Klettfesterr(Sel355P,j)=Btemperr(Sel355P);
            BetaAer355Klettfest(Sel355P,j)=Betaaer(Sel355P); BetaAer355Klettfesterr(Sel355P,j)=Betaaererr(Sel355P,1);
            
            
            % und jetzt rechnen wir noch einmal 532S mit dem Lidarverhaltnis aus LR532P
            % dies gibt die "ausP" Loesungen
            
            
            [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532Sarr(j), FitRangearr(:,j), H(Sel532S), P532SKlett(Sel532S,j), P532SKletterr(Sel532S,j), ...
                LR532arr(Sel532S,j), AlRay532(Sel532S,1), BeRa532S(Sel532S,1));
            Betaaer(Sel532S)=Beta-BeRa532S(Sel532S,1);
            Betaaererr(Sel532S,1)=abs(dBdR.*BSRAtFiterr)+abs(dBdLR.*LR532arrerr(Sel532S,j))+abs(dBdP.*P532SKletterr(Sel532S,j));
            Btemp(Sel532S)=Beta./BeRa532S(Sel532S,1); Btemperr(Sel532S)=Betaaererr(Sel532S,1)./BeRa532S(Sel532S,1);
            BSR532SKlettausP(Sel532S,j)=Btemp(Sel532S); BSR532SKlettausPerr(Sel532S,j)=Btemperr(Sel532S);
            BetaAer532SKlettausP(Sel532S,j)=Betaaer(Sel532S); BetaAer532SKlettausPerr(Sel532S,j)=Betaaererr(Sel532S,1);
            
            
            
        end % for Zeitschritte "entries"
        
        
        
        
        P532roh=P532A';
        P532Sroh=P532SA';
        P355roh=P355A';
        
        speichern532p=['P532Klett P532Kletterr P532roh BSRAtFit532arr LR532arr BetaAer532Klett BetaAer532Kletterr ' ...
            'AlphaAer532 AlphaAer532err BSR532Klett BSR532Kletterr BetaAer532Klettfest BetaAer532Klettfesterr BSR532Klettfest ' ...
            'BSR532Klettfesterr dBeta532dP  dBeta532dR dBeta532dLR Abbruch532 P532Abackground'];
        
        speichern532s=['P532SKlett P532SKletterr P532Sroh BSRAtFit532Sarr LR532Sarr BetaAer532SKlett BetaAer532SKletterr ' ...
            'AlphaAer532S AlphaAer532Serr BSR532SKlett BSR532SKletterr BSR532SKlettfest BSR532SKlettfesterr BetaAer532SKlettfest ' ...
            'BetaAer532SKlettfesterr dBeta532SdP  dBeta532SdR dBeta532SdLR Abbruch532S P532SAbackground BetaAer532SKlettausP BSR532SKlettausP'];
        
        
        speichern355p=['P355Klett P355Kletterr P355roh BSRAtFit355arr  LR355arr BetaAer355Klett BetaAer355Kletterr ' ...
            'AlphaAer355 AlphaAer355err BSR355Klett BSR355Kletterr BSR355Klettfest BSR355Klettfesterr BetaAer355Klettfest ' ...
            'BetaAer355Klettfesterr dBeta355dP  dBeta355dR dBeta355dLR Abbruch355 P355Abackground'];
        
        sonstiges = [' H matlabzeit Wolkenmaske FitRangearr BSR532sollarr BSR355sollarr BSRAtFit532start BSRAtFit355start VolDep532' ];
        
        
        liste=[speichern532p ' ' speichern532s ' ' speichern355p sonstiges];
        
        speichern =['save ' speicherdir Speichernamepraefix filename ' ' liste];
        
        eval(speichern)
        
        
    end % for alle Daten eines Tages
    
end % if files nicht leer



return



