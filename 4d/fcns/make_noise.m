function [out, nstot, svfile, r] = make_noise(params, mode, perc,num_noise, svfile, out, dataclass)
nrecx=params.nrecx; 
nrecy=params.nrecy; 
nsrcx=params.nsrcx;
nsrcy=params.nsrcy;
indm = params.indm;

alpha = params.alpha; 
nstot = zeros(size(out));

if mode.noise==1
    rng(2) 
    svfile = [svfile 'deno']; 
%     tot_noise = 0; 
%     perc = .8; %used to be 0.01
    idx = setdiff(1:nsrcx*nsrcy,indm);%set of guy's you've observed
    indn = randperm(length(idx)); %random permutation of them
    r = floor(length(indn)*perc); %define length for ease
    %here, we corrupt with really strong points of noise  
    for i = 1:r
        switch(dataclass)
            case{'l2'}
%                 alpha = 1e-2;
                normb = norm(out(:,:,idx(indn(i))),2);
                noisemat = alpha*normb*(randn(nrecx, nrecy)+1i*randn(nrecx,nrecy));
            case{'l1', 'l0'}
                %here we want to add very large noise values to a few
                %points
%                 alpha = 1e-1; 
                normb = norm(out(:,:,idx(indn(i))),1); 
                noisemat = vec(alpha*normb*(randn(nrecx, nrecy)+1i*randn(nrecx,nrecy))); %formerly 1e3
                [~, nind] = sort(noisemat, 'descend');
                noisemat(nind(num_noise+1:end))=0; %lose everything but the strongest num_noise
                noisemat=reshape(noisemat, nrecx, nrecy);
            case{'linf'}
%                 alpha = 5e-3; 
                normb = norm(out(:,:,idx(indn(i))),inf); 
                noisemat = alpha*normb*exp(randn([nrecx, nrecy])*1i);
        end
        out(:,:,idx(indn(i))) = out(:,:,idx(indn(i))) + noisemat; %add back to values? 
        nstot(:,:,idx(indn(i))) = noisemat;
    end
%     out(:,:,idx(indn(1:floor(length(indn)*perc))))= out(:,:,idx(indn(1:floor(length(indn)*perc))))...
%         + 1e3*(randn(nrecx,nrecy,length(indn(1:floor(length(indn)*perc))))...
%         +1i*randn(nrecx,nrecy,length(indn(1:floor(length(indn)*perc)))));
else 
%     perc = .8; %used to be 0.01
    idx = setdiff(1:nsrcx*nsrcy,indm);%set of guy's you've observed
    indn = randperm(length(idx)); %random permutation of them
    r = floor(length(indn)*perc); %define length for ease 
    switch(dataclass)
        case{'l2'}
%             alpha = 1e-2;
        case{'l1', 'l0'}
%             alpha = 1e-1; 
        case{'linf'}
%             alpha = 5e-3; 
    end
end

end