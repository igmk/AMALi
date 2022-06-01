function WQ = raybckwq(WaveLength,InPol,OutPol,T,rho)
%RAYBCKWQ  Rayleigh-Rueckstreu-WQ 
%
%            WQ = raybckwq(WaveLength,InPol,OutPol)
%
%          Wavelength : wavelength in m
%          InPol      : ['p','s','u'] 
%                       polarization state of incident light
%          OutPol     : ['p','s','u'] 
%                       detected polarization of scattered light
%
%          Wellenlaenge in m.
%

%	Gebe 27-03-92

  WaveLength = crrctwvl(WaveLength);

  RI = stdriair(WaveLength,T,rho);
  RI2 = RI.^2;
  
  WQ1= pi^2 * ((RI+1).*(RI-1)).^2        ...
        ./ (WaveLength.^2 * stddicht).^2 ...
        .* anisatpi(InPol,OutPol);

 %% diese urspr"ungliche Formel .* 1.1304 um konsistent mit Buchholtz - dann "ahnlich zu 
 %% unterer Formel
 
WQ = 9.*pi^2.*(RI2-1).^2 ./ (WaveLength.^4.*6.4894e50.*(RI2+2).^2) .* anisatpi(InPol, OutPol);
%% nach Mc Cartney, Optics of the Atmosphere Wiley, Neu York, 1976 Kap 4. CR 6/2002

%% wenden hier schon Temperatur-Korrektur an, bislang +15Grad Celsius
%WQ = WQ .* 288.15 ./ T;
WQ = WQ1.*1;     %%%.*0.9;                            %%%%1.1304;