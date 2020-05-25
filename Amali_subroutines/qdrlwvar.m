function Res = qdrlwvar(x,y)
%QDRLWVAR  Bilde Riemann-Summe von x bis MAX(x), x laeuft zwischen
%          MIN(x) und MAX(x).


% 97-02-05 optimized /GeBe


%	(c) GeBe 10-02-92

  myassert(size(y)==size(x))
  myassert(size(x,1)==1 | size(x,2)==1)
  myassert(size(y,1)==1 | size(y,2)==1)

  Res = x * 0;

  Integr       = ( y(1:end-1) + y(2:end)) / 2 .* diff(x);
  Tmp          = cumsum(Integr(end:-1:1));
  Res(1:end-1) = Tmp(end:-1:1);


