function String = fnmerge(Node, Device, Directory, Filename, Extension, Version)
%FNMERGE   assemble complete filename
%
%            String = fnmerge(Node,Dev,Dir,Filename,Ext,Ver)
%

% 14.10.92 gebe

  String = '';

  if nargin < 6,
    Version = [];
  end 

  if nargin < 5,
    Extension = [];
  end 

  if nargin < 4,
    Filename = [];
  end 

  if nargin < 3,
    Directory = [];
  end 

  if strcmp(computer, 'VAX') ...
    | strcmp(computer, 'VAX_VMSG') ...
    | strcmp(computer, 'VAX_VMSD'),

    if length(Node) > 0
      String = [ String Node '::'];
    end

    if length(Device) > 0
      String = [ String Device ':'];
    end

    if length(Directory) > 0
      String = [ String '['];
      for i = length(Directory(:,1)),
        String = [ String deblank(Directory(i,:)) '.'];
      end
      String(length(String)) = ']';
    end

    if length(Filename) > 0
      String = [ String Filename '.'];
    end

    if length(Extension) > 0
      String = [ String Extension];
    end

    if length(Version) > 0
      String = [ String ';' Version];
    end

  elseif opersyst == 'UNIX  ' | opersyst =='???-OS' ,

    if length(Directory) > 0,   
      String = [ String '/'];
      lenD = length(Directory(1,:));  %%(:,1))
      for i = 1:lenD,
        String = [ String char(deblank(Directory(:,i))) '/'];   %%(i,:)
      end
      String(length(String)) = '/';
    end

    if length(Filename) > 0
      String = [ String Filename '.'];
    end

    if length(Extension) > 0
      String = [ String Extension];
    end
    
% *****************************************************
%
%
%
  elseif strcmp(computer,'PCWIN'),

    assert( isempty( Device) | isempty( Node));

    if ~isempty( Device),
      String = [ String Device ':'];
    end

    if ~isempty( Node),
      String = [ String '\\' Node];
    end

    if ~isempty( Directory) > 0
      String = [ String '\'];
      lenD = length(Directory(:,1));
      for i = 1:lenD,
        String = [ String deblank(Directory(i,:)) '\'];
      end
      String(length(String)) = '\';
    else
      String = [ String '\'];    
    end

    if length(Filename) > 0
      String = [ String Filename '.'];
    end

    if length(Extension) > 0
      String = [ String Extension];
    end
    
  else

    error([upper(mfilename) ': OS ' computer 'not yet supported']);

  end


