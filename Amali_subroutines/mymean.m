function [y, yAE] = mymean(x, xAE)
%MYMEAN	Average or mean value.
%	For vectors,  MYMEAN(X,XAE)  is the mean value of the elements in X.
%	For matrices,  MYMEAN(X,XAE) is a row vector containing the mean value
%	of each column.
%
%	See also MEDIAN, STD, MIN, MAX.

%	Copyright (c) 1984-93 by The MathWorks, Inc.
%	Copyright (c) 1998 Georg Beyerle

  [m,n] = size(x);

  Fin = ~isnan( x);
  Idx = find( Fin);

  weight  = str2num(int2str(reshape( Fin, m, n)));
  x0      = weight;
  x0(Idx) = x(Idx);
  if nargout > 1,
    x0AE      = weight;
    x0AE(Idx) = xAE(Idx);
  end
  
  SumWeight = sum( weight);

  Idx            = find( SumWeight==0);
  SumWeight(Idx) = NaN;

  y   = sum( x0) ./ SumWeight;
  if nargout > 1,
    yAE = sqrt( sum( x0AE.^2)) ./ SumWeight;
  end

