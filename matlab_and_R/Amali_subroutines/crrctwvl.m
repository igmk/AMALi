function RealWvl = crrctwvl(Wvl)
%CRRCTWVL correct wavelengths
%
%           RealWvl = crrctwvl(Wvl)
%

% (c) GeBe 1996

if Wvl == 2.82e-7, Wvl = 5.32e-7; end

  MFile = [upper(mfilename) ': '];

  RealWvl = Wvl * NaN;

  WvlList = [ 266.04e-9,  289.06e-9, 299.06e-9, ...
              307.95e-9,  331.77e-9, ...
              337.10e-9, ...
              353.07e-9,  354.71e-9,  384.83e-9, ...
              386.66e-9,  407.49e-9,  523.50e-9, 532.07e-9, ...
              607.35e-9, 660.38e-9, 910.00e-9 1064.14e-9, 443.0e-9, ...
              999e-9];  %%% Default value

  Len = length( WvlList);
  for k = 1:Len,
    Idx = find(abs(Wvl - round(WvlList(k)*1e9)*1e-9) < 1e-9 | ...
               abs(Wvl - WvlList(k)) < 2e-11);
    myassert( isnan( RealWvl(Idx)));
    RealWvl(Idx) = WvlList(k);
  end
  
  NaNWvl = isnan(RealWvl);
  if any(NaNWvl),
    Idx = find(NaNWvl);
    disp([upper(mfilename) ': dont know how to correct wavelength ' ...
           num2str(Wvl(Idx(1))*1e9) ' nm'])
%    warning('***')
    RealWvl(Idx) = Wvl(Idx);
  end
  