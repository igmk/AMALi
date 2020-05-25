function C = concat(Dim,A,B)
%CONCAT concatenate two matrices or cell arrays of different dimensions
%
%         C = concat(Dim,A,B)
%
%         Dim = 1 : concatenate row wise
%         Dim = 2 : concatenate column wise
%

% (c) GeBe 1996

  assert( nargin == 3)
  assert( iscell(A) == iscell(B))

  [RowA,ColA] = size(A);
  [RowB,ColB] = size(B);

  if iscell( A) | isstruct( A),
	  if Dim == 2,
	    C = [A B];
	  elseif Dim == 1,
	    C = [A; B];
	  else
	    error([upper(mfilename) ': not yet implemented']);
	  end
  else
	  if Dim == 2,
	    C = ones( max( RowA, RowB), ColA+ColB) * NaN;
	    C( 1:RowA, 1:ColA)        = A;
	    C( 1:RowB, ColA+(1:ColB)) = B;
	  elseif Dim == 1,
	    C = ones( RowA+RowB, max( ColA, ColB)) * NaN;
	    C( 1:RowA, 1:ColA)        = A;
	    C( RowA+(1:RowB), 1:ColB) = B;
	  else
	    error([upper(mfilename) ': not yet implemented']);
	  end
  end  

