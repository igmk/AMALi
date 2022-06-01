function WQ = raytotwq(WaveLength, T, rho)
%RAYTOTWQ  Rayleigh-WQ als Funktion der Wellenlaenge.
%
%             WQ = raytotwq(WaveLength)
%
%          Wellenlaenge in m.

%	Gebe 27-03-92

  WaveLength = crrctwvl(WaveLength);

  RI = stdriair(WaveLength, T, rho);
  RI2 = RI.^2;
  
  WQ1 = 8 * pi / 3 ...
        * pi^2 * ((RI+1).*(RI-1)).^2 ...
        ./ (WaveLength.^2 * stddicht).^2 ...
        * aniso4pi;



WQ = 24.*pi^3.*(RI2-1).^2 ./ (WaveLength.^4.*6.4894e50.*(RI2+2).^2) .*aniso4pi;
WQ=WQ1.*1;     %%.*0.9;   Tests    
%% korrigieren hier f"ur Temp. (bislang +15C
%WQ = WQ .* 288.15./T;
%WQ = 1.1304.*WQ;
 
 