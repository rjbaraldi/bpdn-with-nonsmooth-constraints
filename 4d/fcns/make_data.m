function [out, indm, svfile] = make_data(params, mode, perc, out, svfile)

if mode.interp==1
    rng(1)
    index = randperm(params.nsrcx*params.nsrcy);%generate random index for 
    indm  = index(1:floor(params.nsrcx*params.nsrcy*perc));%take the first
% loind = indm(1); 
% left_out = squeeze(out(:,:,loind)); 
    out(:,:,indm) = 0;
    svfile = [svfile, 'inter']; 
else
    indm = [] ;
    
end



end