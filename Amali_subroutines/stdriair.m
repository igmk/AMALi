function Res = stdriair(WaveLength, T, rho)
%STDRIAIR Brechungsindex von Luft bei Standard-Bedingungen.
%         Wellenlaenge in m
%
%           Res = stdriair(WaveLength)

%	Gebe 27-03-92

  WLSqr = WaveLength.^2;

  Res = ( 6432.8 ...
           + ( 2949810.0) ./ (146.0 - 1.0e-12 ./ WLSqr) ...
	   + ( 25540.0)   ./ ( 41.0 - 1.0e-12 ./ WLSqr)) ...
	 * 1.0e-8;

%
% Umrechnen auf Std-Temperatur
%
  Takt = stdtempe-20;
  Res = Res * (stdtempe + 15) / Takt;
  %Res = Res.*1.4;
  Res = Res + 1;

 % return;
  
% ******************************************************

%
%  Anthony Buchholz  Applied Optics 1995, Vol 34; No 15; S. 2765 ff
  

ll = WaveLength .* 1e6;  %% in Mikrometern
T1 = 5791817./(238.0185 - (1./ll).^2) + 167909 ./ (57.362 - (1./ll).^2);
Res = T1./1e8; %% gilt f"ur 15Grad C und 1013,25hPa
%Res = Res.*1.0.*rho./rho(1);
Res = Res +1;
  

return

% ********************************************************

%
%  Parametrization by Edlin (1950 ?)
%
  l0 = WaveLength * 1e6;
  N = 3.*(64.328 + 29498.1 ./ (146 - 1 ./ l0.^2) + 255.4 ./ (41 - 1 ./ l0.^2))*1e-6;

  H2OPartialPress = 0;
  H2OPartialPress = 2300;
  Temperature = 300;

  t = Temperature - stdtempe;
  t = stdtempe-15;   %%%T;
  f = H2OPartialPress * 760 / 101325;  % mmHg -> Pa
  N = N - (0.0624 - 0.00068 / l0^2) ./ (1 + 0.003661 * t) * f * 1e-6;

  
  %%%N = N * (stdtempe + 15) / stdtempe;
  N = N * (stdtempe) /(stdtempe-15);
  
  N = N + 1;
  
  Res = N;
  
  return;
  
% **************************************************  
%
%  Parametrization by Ciddor, AO 35, p1566 (1996)
%
  sig = 1 ./ WaveLength * 1e-6;

  k0 = 238.0185;
  k1 = 5792105;
  k2 = 57.362;
  k3 = 167917;

  Nas = 1 + (k1 ./ (k0 - sig.^2) + k3 ./ (k2 - sig.^2))*1e-8;
  
  (Nas^2 - 1)^2
  
  xc = 450;
  Naxs = 1 + (Nas-1) * (1+0.534e-6*(xc-450));
  
  (Naxs^2 - 1)^2

u = 1;


% *******************************************************
%
% The Refractive Index of Air B Edlen 1966
% http://iopscience.iop.org/0026-1394/2/2/002
%
% s Wellenzahl in Mikrometern
% p Druck in Torr 133,322 Pa

s = 1./(WaveLength.*1e6);
t = T +273.15;
p = rho.*T.*boltzman./133.322;

n1h = 8342.13 + 2406030.*(130-s.^2).^(-1) + 15997.*(38.9-s.^2).^(-1);
n1tp = n1h.*0.00138823 .*p./(1 + 0.003671 .*t);

Res = n1tp./1e8+1;
%%% das gibt so viel zu kleine AlRay ....
Res = Res(1,1);



% *******************************************************
%
% E. D. Peck and K. Reeder, J. Opt. Soc. Am. 62, 958 (1972)
%
%
% stimmt fuer UV mit Edlen 1966 ueberein  

n1h = 5791817 ./(238.0185 -s.^2) + 167909 ./(57.362-s.^2);
 
n1tp = n1h.*0.00138823 .*p./(1 + 0.003671 .*t);
Res = n1tp./1e8+1;




