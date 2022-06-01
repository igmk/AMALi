function ConnRangeIdx = connrnge(X,Sort)
%CONNRNGE  range of sequences of 1's in vector or matrix
%
%            ConnRangeIdx = connrnge(X,Sort)
%
%          calculate start and end index of connected ranges
%          and sort according to size
%
%            Sort = 1: sort according to size
%            Sort = 0: dont sort
%

% (c) 1996 GeBe

  if nargin < 2,
    Sort = 0;
  end

%
%  only zeros or ones allowed
%
  myassert(X == 1 | X == 0);
  myassert(~any(isnan(X)))

%
%  if row vector, change to column vector
%
  [r,c] = size(X);

  if r == 1,
    X = X(:);
    [r,c] = size(X);
  end 

%
%  add zeros to start and end
%
  X = [zeros(1,c);X;zeros(1,c)];
  ConnRangeIdx = [];

  for i = 1:c,

    Df = diff(X(:,i) ~= 0);

    Start = find(Df == 1);
    End   = find(Df == -1) - 1;

    if Sort,

      [Tmp,SortIdx] = sort(End - Start);
      SortIdx       = SortIdx(length(SortIdx):-1:1);
      Tmp           = [Start(SortIdx)'; End(SortIdx)'];

    else

      Tmp = [Start'; End'];

    end

    Tmp          = Tmp(:);
    ConnRangeIdx = concat(2,ConnRangeIdx,Tmp);

  end


