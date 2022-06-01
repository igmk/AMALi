function res = anisatpi(InPol,OutPol)
%ANISATPI  Anisotropiefaktor bei Streuung in Rueckwaertsrichtung.
%          Aufruf
%
%            res = anisatpi(InPol,OutPol)
%

%	Gebe 27-03-92

  myassert(size(InPol) == size(OutPol))

  Len = length(InPol);
  
  res = ones(size(InPol)) * NaN;

  for i = 1:Len

    if     InPol(i) == 'u' & OutPol(i) == 'u',
%    res =  0.5 * (( 48 * epsiluft + ( 180 + 4 * epsiluft) * 2 ) / 180);
      res(i) = ( 180 + 28 * epsiluft ) / 180;   %%%%% cr 12/2009 hoppala sehe hier eine 90 im Nenner
    elseif InPol(i) == 'p' & OutPol(i) == 'u',
      res(i) = ( 180 + 28 * epsiluft ) / 180;
    elseif InPol(i) == 'p' & OutPol(i) == 's',
      res(i) = ( 12. * epsiluft ) / 180;
    elseif InPol(i) == 'p' & OutPol(i) == 'p',
      res(i) = ( 180 + 16 * epsiluft ) / 180;
    else
      error([upper(mfilename) ': ungueltige Kombination von Polarisationen : ' ...
            InPol ' ' OutPol]);
    end
  end

