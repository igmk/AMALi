function Res = qdrupvar(x,y)
%QDRLWVAR  Bilde Riemann-Summe von x bis MAX(x), x laeuft zwischen
%          MIN(x) und MAX(x).
%% Achtung Verwendung Res=qdrupvar(dx,f(x))

% 97-02-05 optimized /GeBe


%	(c) GeBe 10-02-92
% cr 2/2013   behandle NaN im Integral als 0

if size(x) ~= size(y),
    if length(x) == length(y), y=y'; end
end

  Res = x * 0;

 if size(x) ~= size(y)
     return
 end
  

 if length(x) >2
  Integr     = ( y(1:end-1) + y(2:end)) / 2 .* diff(x);
  wo = find(isnan(Integr)); Integr(wo)=0;
  Res(2:end) = cumsum(Integr);
 end

