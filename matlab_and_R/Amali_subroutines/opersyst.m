function OS = opersyst()

  MFile = 'OPERSYST: ';

  if strcmp( computer, 'VAX') ...
     | strcmp( computer, 'VAX_VMSG') ...
     | strcmp( computer, 'VAX_VMSD'),
    OS = 'VAXVMS';
  elseif  strcmp( computer, 'PC') ...
        | strcmp( computer, '386') ...
        | strcmp( computer, 'PC386'),
    OS = 'MS-DOS';
  elseif  strcmp( computer, 'PCWIN'),
    OS = 'MS-WIN';
  elseif  strcmp( computer, 'SUN4') ...
        | strcmp( computer, 'sparc-sun-sunos4.1.3_U1') ...
        | strcmp( computer, 'sparc-sun-sunos4.1.3'),
    OS = 'UNIX  ';
  else
    OS = '???-OS';
  end

  assert(length(OS)==6)

