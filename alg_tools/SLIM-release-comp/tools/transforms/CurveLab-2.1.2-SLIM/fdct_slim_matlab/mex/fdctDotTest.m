p = path;path(p,'/Users/crussell/CurveLab-2.0-ext/fdct_wrapping_cpp/mex/');
p = path;path(p,'/Users/crussell/CurveLab-2.0-ext/nfdctm_wrapping_cpp/mex/');

m = 512;
n = 500;
CLini(m,n,1,6,16);
  global HDR
  global CN
  global M
  global N

  global AC
  global NBSCALES
  global NBANGLES_COARSE

x = randn(m,n);
Cx = CL(x);
y = rand(size(Cx,1),1);
CHy = CL(y,'transp');

disp('<x,CHy> = ');
A = norm( sum(x(:).*conj(CHy(:))) )
disp('<Cx,y> = ');
B = norm( sum(Cx(:).*conj(y(:))) )
