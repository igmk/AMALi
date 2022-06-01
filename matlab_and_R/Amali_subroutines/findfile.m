function Files = findfile(SuchString)
%FINDFILE   search directory for files
%
%             Files = findfile(SuchString)
%
%           Liefert Textmatrix aller Files in einem Verzeichnis,
%           die Suchmuster 'SuchString' entsprechen.
%           Liefert [], falls keine Match.
%

% (c) GeBe 1997

  assert( size(SuchString,1) == 1);
  
  List = dir(SuchString);

%
%  Bug in 'dir' command 
%
  if findstr(SuchString,'*')
    [Node,Dev,Dir,Name,Ext,Ver] = fnsplit(SuchString);
    Directory = fnmerge(Node,Dev,Dir,[],[],[]);
  else
    Directory = [];  
  end
  
  Files = [repmat(Directory,length(List),1) strvcat(List.name)];

%
%  'dir' doesn't sort
%
  Files = sortrows(Files);

