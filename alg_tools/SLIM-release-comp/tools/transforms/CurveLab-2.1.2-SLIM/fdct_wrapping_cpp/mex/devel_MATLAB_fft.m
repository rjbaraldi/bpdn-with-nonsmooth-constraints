function devel_MATLAB_fft()
    Sm=32; Pm=Sm/2; Pm1=Pm; Pm2=Pm+2;
    Sn=64; Pn=Sn/2; Pn1=Pn; Pn2=Pn+2;
    Sf=sqrt(Sm*Sn);
    nbscales = floor(log2(min(Sm,Sn)))-3;
    nbangles_coarse = 16;
    allcurvelets = 0;
    fprintf('size=(%d,%d) nbscales=%d nbangles_coarse=%d allcurvelets=%d\n',Sm,Sn,nbscales,nbangles_coarse,allcurvelets)
    fprintf('Pm=%d:%d Pn=%d:%d\n',Pm1,Pm2,Pn1,Pn2)
  
    % make input
    X=rand(Sm,Sn);
    x=fft2(X);
    fx=x(1:3,1:3)'
    sx=fftshift(x)/Sf;
    sfx=sx(Pm1:Pm2,Pn1:Pn2)'

    % opCurvelet
    opC = opCurvelet(Sm,Sn,nbscales,nbangles_coarse,allcurvelets);
    oC = opC*X(:);
    oY = opC'*oC;
    opn=norm(X(:)-oY(:));

    %call full function
    C = fdct_wrapping_mex(Sm,Sn,nbscales, nbangles_coarse, allcurvelets, double(X));
    Y = ifdct_wrapping_mex(Sm,Sn,nbscales, nbangles_coarse, allcurvelets, C);
    on=norm(X(:)-Y(:));
    oN=norm(oY(:)-Y(:));

    % opCurveletNoFFT
    opCN = opCurveletNoFFT(Sm,Sn,nbscales,nbangles_coarse,allcurvelets);
    oCN = opCN*sx(:);
    osy = opCN'*oCN;
    oCNn=norm(oC-oCN);
    opnN=norm(sx(:)-osy(:));
    rsfy=reshape(osy,size(sx));
    osfy=rsfy(Pm1:Pm2,Pn1:Pn2)'
    osy=Sf*ifftshift(rsfy);
    ofy=osy(1:3,1:3)'

    %call nofft function
    c = fdct_wrapping_nofft_mex(Sm,Sn,nbscales, nbangles_coarse, allcurvelets, double(sx));
    sy = ifdct_wrapping_nofft_mex(Sm,Sn,nbscales, nbangles_coarse, allcurvelets, c);
    sn=norm(sx(:)-sy(:));
    sfy=sy(Pm1:Pm2,Pn1:Pn2)'
    sy=Sf*ifftshift(sy);
    fy=sy(1:3,1:3)'
    y=ifft2(sy);
    fn=norm(X(:)-y(:));
    fN=norm(oY(:)-y(:));

    fprintf('norms:\topn=%f orig=%f/%f\n\topNn=%f/%f fft=%f fin=%f/%f\n',opn,on,oN,oCNn,opnN,sn,fn,fN)
end
