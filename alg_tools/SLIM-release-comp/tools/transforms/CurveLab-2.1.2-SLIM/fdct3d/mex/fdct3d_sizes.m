function C = fdct3d_sizes(m,n,p)
  nbscales = floor(log2(min([m,n,p])))-2;
  nbdstz_coarse = 8;
  allcurvelets = 0;
  C = fdct3d_sizes_mex(m,n,p,nbscales,nbdstz_coarse,allcurvelets);