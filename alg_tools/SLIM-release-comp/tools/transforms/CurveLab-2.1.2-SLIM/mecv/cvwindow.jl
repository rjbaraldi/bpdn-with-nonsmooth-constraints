function [wl,wr] = cvwindow(x)
  wr = zeros(Float64,size(x));
  wl = zeros(Float64,size(x));
  
  eps = 1e-16;
  sml = find(y->y<eps,x); %too small
  mid = find(y->y>=eps & y<=1-eps,x); %just right
  lrg = find(y->y>1-eps,x); %too large
  
  wl[sml] = 0;  wr[sml] = 1;
  wl[lrg] = 1;  wr[lrg] = 0;
  
  xmid    = x[mid];
  a       = exp.(1-1./(1-exp.(1-1./(1-xmid))));
  b       = exp.(1-1./(1-exp.(1-1./xmid)));
  n       = sqrt(a.^2 + b.^2);
  wl[mid] = a./n;  wr[mid] = b./n;
  return wl, wr
end