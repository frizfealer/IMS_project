function [ X ] = RandInplacePermute( X )
%RandInplacePermute A faster version of randpermutation
%this is Fisher¡VYates shuffle algorithm
for i = 1:numel(X)
  w = ceil(rand * i);
  t = X(w);
  X(w) = X(i);
  X(i) = t;
end

end

