function nq=nquist(n)

%% n=even:
% returns the index of +nyquist for n even
% n=odd: nyquist is not represented, it is half way between (n+1)/2
% and (n+1)/2+1, so the sample relating to the highest positive frequency 
% is returned (=(n+1)/2).

    if mod(n,2) ~= 0
        nq=(n+1)/2;
    else
        nq=n/2+1;
    end

end
