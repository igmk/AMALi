function Res = strstr(StringMatrix,Pattern)
%STRSTR	Textmatrix nach String absuchen.  
%       STRSTR(StringMatrix,Pattern) sucht jede Zeile der Stringmatrix 
%       'StringMatrix' nach dem ersten Vorkommen von 'Pattern' ab, und 
%       liefert bei Erfolg den Index in der Zeile im zurueck. Taucht
%       'Pattern' in der Zeile nicht auf, ist das Funktionsergebnis 0.

%	GeBe 21-04-92
%            03-09-94 : Tuning /gebe
%	Copyright (c) 1992 by GeBe

  MFile = 'STRSTR: ';

  if nargin ~= 2,
    error([MFile 'Aufruf ist Res = strstr(String,Pattern)']);
  end

%
% Konvertiere Strings in Integers und invertiere StringMatrix
%
  Pattern      = abs(Pattern);  
  StringMatrix = abs(StringMatrix');  

  [RowStr ColStr] = size(StringMatrix);
  [RowPat ColPat] = size(Pattern);

  SMLen = RowStr * ColStr;

  if RowPat ~= 1,
    error([MFile 'als Pattern sind nur Strings zulaessig!']);
  end

  if ColPat > RowStr,
    Res = zeros(ColStr,1);
    return;
  end

%
% um  'xxxxPat'  zu verhindern.
%     'ternxxx'
%
  if ColPat > 1,
    t = [zeros(ColStr,ColPat-1) ones(ColStr,RowStr-ColPat+1)]';
    t = t(:);
  else
    t = ones(SMLen,1);
  end

  StringMatrix = StringMatrix(:);

  v = StringMatrix(ColPat:SMLen) == Pattern(ColPat) & t(ColPat:SMLen);

  for i = ColPat-1:-1:1,
    v = StringMatrix(i:SMLen-ColPat+i) == Pattern(i) & v;
  end

%
%  Forme um:
%
%   [ 0 1 0 0 1 0 0              [ 0 6 0 0 3 0 0              
%     0 0 0 0 0 0 0                0 0 0 0 0 0 0
%     1 0 0 1 0 1 0      ->        7 0 0 4 0 2 0
%     0 1 0 0 0 1 0 ]              0 6 0 0 0 2 0 ]
%

  if ColPat > 1,
    v = [v; zeros(ColPat-1,1)];
  end

  v = reshape(v,RowStr,ColStr) .* ((RowStr:-1:1)' * ones(1,ColStr));

%
%  liefert mit obiger Matrix :  [6 0 7 6],[2 1 1 2]
%

  [Max, MaxIdx] = max(v);

  Res = ((~(~Max)) .* MaxIdx)';

