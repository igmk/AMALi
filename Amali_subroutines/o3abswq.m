function [WQ,WQAE] = o3abswq(WaveLength,Temperatur)
%O3ABSWQ Ozon-Absorptions-WQ als Funktion der Wellenlaenge.
%        Wellenlaenge in m.
%
%          [WQ,WQAE] = o3abswq(WaveLength,Temperatur)
%

%
% Note:
%
%  Obtain x section from Bass/Paur. Use quadratic fit given by
%  Bass/Paur, as their fit is a weighted fit and includes varying
%  accurracies of the data points. 
%

%
% Note:
%
%  The DIAL wavelengths are: 307.95 and 353.15 nm
%

%	Gebe 27-03-92
%  einige konkrete Zahlen gefunden, die auf Brasseur Solomon 1984 zurueck
%  gehen  cr 2/2015
% beste Quelle sei Landolt-Boernstein 1988 Springer Band 4b physical and
% chemical properties of the air, pages 526-529
% Nicht implementiert, aber Brion et al. 1998 ist auch eine Quelle
% http://igaco-o3.fmi.fi/ACSO/files/brion_et_al_1998.pdf
% dort ca 2.48-21 (cm.^2) bei 532nm, Nachweis das Temperaturabhaenigkeit
% nur schwach


  WaveLength = crrctwvl(WaveLength);

  WQ   = 0;
  WQRE = 0.05;   % ???

%
%  TMF - Routine:
%
%               IF (WAVELENGTH.EQ.307.95) THEN  
%                  CROSS_SECTION_OUT(J,I) = 
%     *               (1.2853E-4)*TCELSIUS*TCELSIUS+
%     *               (2.7569E-2)*TCELSIUS+
%     *               (12.857)
%              ELSEIF (WAVELENGTH.EQ.353.15) THEN
%                  CROSS_SECTION_OUT(J,I) = 1.0E-3
%               ELSEIF (WAVELENGTH.EQ.331.77) THEN
%                  CROSS_SECTION_OUT(J,I) = 
%     *               (1.2744E-5)*TCELSIUS*TCELSIUS+
%     *               (3.1243E-3)*TCELSIUS+
%    *               (4.3592E-1)
%               ELSEIF (WAVELENGTH.EQ.384.83) THEN
%                 CROSS_SECTION_OUT(J,I) = 1.0E-3
%              ELSE
%                  STOP 'CANNOT CALCULATE CROSS SECTION'
%               ENDIF  
%               CROSS_SECTION_OUT(J,I) = CROSS_SECTION_OUT(J,I)*1.0E-20

  if nargin < 2,
    Temperatur = stdtempe;
  end

%  DataTemp = [-70 -55 -45 -30 0 25] + stdtempe;   % use Kelvin
  DataTemp = [-70 -55 -45 -30 0 25];   % use Celsius

  Data1 = [ ...
   283.609   275.45 279.60 279.37 283.44 287.45 288.57;
   283.659   272.15 276.15 275.72 278.35 283.19 285.95;
   285.859   213.96 217.93 218.87 218.30 214.82 225.51;
   285.909   213.42 216.91 216.77 218.77 213.89 223.44;
   288.959   146.90 149.64 150.78 152.27 152.67 157.73;
   289.009   146.77 149.28 149.54 151.36 152.31 158.03;
   289.859   132.46 135.27 135.03 137.29 138.40 142.16;
   289.909   132.31 134.72 134.58  NaN   137.36 141.85;
   294.158    77.94  78.90  79.50  80.44  81.73  84.65;
   294.208    77.49  78.47  78.95  79.12  81.35  84.60;
   299.052    40.48  40.65   NaN   41.73  43.29  44.80;   
   299.102    39.24  39.96   NaN   41.03  42.40  44.61;   
   307.950    11.49  11.78  11.98  12.08  12.80  13.66;
   353.150     1e-3   1e-3   1e-3   1e-3   1e-3   1e-3;
  ];

  Data2 = [ ...
   283.609  286.86  9.3685e-2  -9.0636e-4;
   283.659  282.89  1.2955e-1  -2.5539e-4;
   285.859  220.74  1.4092e-1   8.4542e-4;
   285.909  219.54  1.1572e-1   5.7993e-4;
   288.959  154.81  1.0722e-1  -6.6120e-6;
   289.009  154.21  1.3344e-1   4.6588e-4;
   289.859  139.64  9.6251e-2  -1.8182e-5;
   289.909  138.54  1.1571e-1   4.5894e-4;
   294.158    NaN      NaN         NaN;
   294.208    NaN      NaN         NaN;
   299.052   43.21  5.9796e-2   2.9783e-4;   
   299.102   42.64  7.0696e-2   3.3370e-4;   
   307.950   12.88  2.7610e-2   1.1424e-4;
   353.150    1e-3  0           0;
  ];

 

  DataWvl = Data1(:,1)   * 1e-9;
  DataXO3 = Data1(:,2:7);
  DataPP  = Data2(:,4:-1:2);

  [Tmp Col] = min(abs(DataWvl - WaveLength));

  if abs(DataWvl(Col)-WaveLength) < 0.6e-9,

    XSect = DataXO3(Col,:);
    pp    = DataPP(Col,:);

    Idx = find(isnan(XSect));
    XSect(Idx) = [];
    DataTemp(Idx)  = [];

    if all(isnan(pp)),
      warning([upper(mfilename) ': calculate fit!'])
      pp = polyfit(DataTemp,XSect,2);
    end
  
    WQ = polyval(pp,Temperatur-stdtempe) * 1e-24;
    WQAE = WQ .* WQRE;
    
    %% neu gesch"atzt f"ur 355nm gem"a/3 obiger Formel
 elseif abs(WaveLength - 355e-9) < 1e-9,
     
     %WQ =1e-25;  %%Rest Hartley Kontinuum cr 6/2005
     WQ = 1.05e-26; % nach Brassuer Solomon 1984 (gefunden in Sica 2001 Ozone corrections for Rayl-Scatter)
     %cr 2/2015


  elseif abs(WaveLength - 607e-9) < 1e-9,

    WQ = 5.2e-25;

  elseif abs(WaveLength - 532e-9) < 1e-9,

    %WQ = 2.81e-25;
    %woher kommt diese Zahl?
    WQ = 2.2e-25; % nach Brassuer Solomon 1984 (gefunden in Sica 2001 Ozone corrections for Rayl-Scatter)

  elseif abs(WaveLength - 523.5e-9) < 1e-9,

    WQ = 2.81e-25;
 
  elseif abs( WaveLength - 110e-9) < 1e-9 | ...
         abs( WaveLength - 385e-9) < 1e-9 | ...
         abs( WaveLength - 387e-9) < 1e-9 | ...
         abs( WaveLength - 380e-9) < 1e-9 | ...
         abs( WaveLength - 390e-9) < 1e-9 | ...
         abs( WaveLength - 391e-9) < 1e-9 | ... 
         abs( WaveLength - 407e-9) < 1e-9 | ...
         abs( WaveLength - 417e-9) < 1e-9 | ...
         abs( WaveLength - 855e-9) < 1e-9 | ...
         abs( WaveLength - 910e-9) < 1e-9 | ...
         abs( WaveLength - 1064e-9) < 1e-9,

    WQ = 1e-28;

  elseif abs(WaveLength - 332e-9) < 1e-9,

    TCELSIUS = Temperatur + stdtempe;
    WQ = 1.2744E-5*TCELSIUS.^2 + 3.1243E-3*TCELSIUS + 4.3592E-1;
    WQ = WQ * 1e-24;

  else

    error([upper(mfilename) ': wavelength ' sprintf('%.2f',WaveLength*1e9) ' nm not implemented'])

  end


