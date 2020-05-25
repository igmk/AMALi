function [Node, Device, Directory, Filename, Extension, Version] = fnsplit(String)
%FNSPLIT   trenne kompleten Pfadnamen
%        
%            [Node,Dev,Dir,Name,Ext,Ver] = fnsplit(String)
%

% 14.10.92 gebe

  MFile = [upper(mfilename) ': '];

  Node      = [];
  Device    = [];
  Directory = [];
  Filename  = [];
  Extension = [];
  Version   = [];

  assert(length(String(:,1)) == 1);

%
%  VMS
%
  if strcmp(computer, 'VAX') ...
    | strcmp(computer, 'VAX_VMSG') ...
    | strcmp(computer, 'VAX_VMSD'),
%
% Node
%
    Idx = strstr(String,'::');
    if Idx > 0
      Node = String(1:Idx-1);
      String(1:Idx+1) = [];
    else
      Node = [];
    end

%
% Device
%
    Idx = strstr(String,':');
    if Idx > 0
      Device = String(1:Idx-1);
      String(1:Idx) = [];
    else
      Device = [];
    end

%
% Directory
%
    if String(1) == '[',
      Idx = strrstr(String,']');
      assert(Idx > 0);
      TmpStr = String(2:Idx-1);
      TmpStr(find( TmpStr == ']')) = [];
      TmpStr(find( TmpStr == '[')) = [];
      Directory = mystrtok(TmpStr,'.');
      String(1:Idx) = [];
    else
      Directory = [];
    end

%
% Version
%
    Idx = strstr(String,';');
    if Idx > 0
      Version = String(Idx+1:length(String));
      String(Idx:length(String)) = [];
    else
      Version = [];
    end

%
% Extension
%
    Idx = strstr(String,'.');
    if Idx > 0
      Extension = String(Idx+1:length(String));
      String(Idx:length(String)) = [];
    else
      Extension = [];
    end

%
% FileName
%
    if length(String) > 0,
      Filename = String;
    else
      Filename = [];
    end

    
       %%%Unix strcmp
       
elseif opersyst == 'UNIX  ' | opersyst == '???-OS',

%
% Node
%
    Idx = strstr(String,'::');
    if Idx > 0
      Node = String(1:Idx-1);
      String(1:Idx+1) = [];
    else
      Node = [];
    end

%
% Device
%
    Idx = strstr(String,':');
    if Idx > 0
      Device = String(1:Idx-1);
      String(1:Idx) = [];
    else
      Device = [];
    end

%
% Directory
%
    if String(1) == '/',
      Idx = strrstr(String,'/');
      assert(Idx > 0);
      TmpStr = String(2:Idx-1);
      TmpStr(find( TmpStr == ']')) = [];
      TmpStr(find( TmpStr == '[')) = [];
      Directory = mystrtok(TmpStr,'/');
      String(1:Idx) = [];
    else
      Directory = [];
    end

%
% Version
%
    Idx = strstr(String,';');
    if Idx > 0
      Version = String(Idx+1:length(String));
      String(Idx:length(String)) = [];
    else
      Version = [];
    end

%
% Extension
%
    Idx = strstr(String,'.');
    if Idx > 0
      Extension = String(Idx+1:length(String));
      String(Idx:length(String)) = [];
    else
      Extension = [];
    end

%
% FileName
%
    if length(String) > 0,
      Filename = String;
    else
      Filename = [];
    end


% ************************************
%
%   MS-Windows/MS-DOS
%
  elseif strcmp(computer,'PCWIN'),

%
% Device
%
    if String(2) == ':',
      Device = String(1);
      String(1:2) = [];
    end

    List = mystrtok(String,'\');
    
    if length( String) >= 2 & strcmp( String(1:2), '\\'),
      Node    = List{1};
      List(1) = [];
    end
    
%
% Directory
%
    Idx = findstr( String, '\');
    if length( Idx) > 1,
      if String(end) == '\',
        Directory = strvcat(List{:});
      else
        Directory = strvcat(List{1:end-1});
      end
    end

%
% Extension and Filename
%
    if String( end) ~= '\',
      Idx = findstr( List{end}, '.');
      if length( Idx) > 0,
        Idx = Idx(end);
        Filename  = List{end}(1:Idx-1);   
        Extension = List{end}(Idx+1:end);   
      else
        Filename = List{end};
      end
    end
   
  else
    error([MFile 'OS ' computer ' not yet supported']);
  end

