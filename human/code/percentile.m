function v = percentile(d,n)


[nr,nc] = size(d);

if nr == 1
  d = d';
  nr = nc;
  nc = 1;
end

x = sort(d);
v = x(1+floor(n*nr),:);
