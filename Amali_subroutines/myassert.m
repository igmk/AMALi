function myassert( Flag, Text)
%ASSERT  teste Bedingung

% gebe 6.7.92

  MFile = 'ASSERT: ';

  if nargin ~= 1 & nargin ~= 2,
    error([MFile 'Aufruf ist assert( Flag[, Text])']);
  end

%  if ~(all(all(Flag(~isnan(Flag))))),
%    if nargin == 1,
%      error([MFile 'assertion failed!']);
%    else
%      error([MFile 'assertion failed: ' Text]);
%    end
%  end

  if ~(all(all(Flag(~isnan(Flag))))),
    dbstop if error
    if nargin == 1,
       error( [MFile 'assertion failed!']);
    else
      error( [MFile 'assertion failed: ' Text]);
    end    
  end

