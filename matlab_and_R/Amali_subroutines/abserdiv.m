function cAE = abserdiv(a,aAE,b,bAE)
%ABSERDIV absolute error of a/b
%
%           cAE = abserdiv(a,aAE,b,bAE)
%
  
% (c) GeBe 1996

  Idx = find(b == 0);
  b(Idx) = NaN * Idx;

%
%  sqrt((aAE/a)^2+(bAE/b)^2) * (a/b)
%
  if length(aAE) > 0 & length(bAE) > 0,
    cAE = sqrt((aAE./b).^2 + (bAE./b.*a./b).^2);
  else
    cAE = [];
  end
