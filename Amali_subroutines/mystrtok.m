function Res = mystrtok(String,Pattern)
%MYSTRTOK String zerlegen
%
%           mystrtok(String,Pattern) 
%
%         sucht den String 'String' nach den Trennzeichen 'Pattern' ab, 
%         und erzeugt einen Text cell array, die die Substrings enthaelt, 
%         zurueck.

%	GeBe 21-04-92
%	(c) 1992 by GeBe
%
%	  97-01-17 : use v5.0 features /GeBe 

  [RowStr ColStr] = size(String);
  [RowPat ColPat] = size(Pattern);

  if RowStr ~= 1 | RowPat ~= 1,
    error([upper(mfilename) ': Argumente keine Strings!']);
  end

  for i=1:ColPat,
    Idx         = find(String == Pattern(i));
    String(Idx) = 0; 
  end

  Idx = connrnge(String > 0, 0);
  Len = length(Idx);

  if Len < 1,
    Res = [];
    return
  end

  k = 1;
  for i = 1:2:Len,
    Res(k) = {String(Idx(i):Idx(i+1))};
    k = k + 1;
  end

