function [maxk, ind]=maxk2(A, k, a_or_b_option)
%note that this produces logical ind
% checkkvalid = (k<=numel(A))&(k>0); 
%returns maximum k elements
if ~isreal(A)
    Aabs = abs(A); 
else
    Aabs = A; 
end

%get the kth largest element
kth = select(Aabs, 1, numel(A), k);

maxk = A(Aabs>=kth);  
if strcmp(a_or_b_option, 'a')
    ind = (Aabs>=kth); 
else
    ind = (Aabs<kth); 
end



end






function kth_ele = select(A, left, right, k)
%returns the k-th smallest element of the list within left/right inclusive
% ie (left <=k<= right)
%the search space within the array is changing for each round, but the list
%is still the same size - therefore k does not need to be updated with each
%round
%slightly faster but fixed memory
while 1
    if left == right
        kth_ele = A(k);
        break
    else
        pindex = median([left, floor((right-left)/2), right]);
%         pindex = medofmed(A, left, right); 
%         pindex = randi(numel(A)); 
        [pindex,A] = partition(A, left, right, pindex); 
        if k==pindex
            kth_ele = A(k); 
            break
        elseif k<pindex
            right = pindex -1;
        else
            left = pindex+1; 
        end
        
    end
    
end
%slightly slower but possibly less memory
% if left==right
%     kth_ele = A(left); 
% else
% %pick a random pivot index
%     pindex = median([left, floor((right-left)/2), right]); 
% %     pindex = medofmed(A, left, right); 
%     [pindex,A] = partition(A, left, right, pindex); 
%     if k==pindex
%         kth_ele = A(k); 
%     elseif k<pindex
%         kth_ele = select(A, left, pindex-1, k);
%     else
%         kth_ele = select(A, pindex+1, right,k); 
%         
%        
%     end
% 
% end
end


function [storeindex, A] = partition(A, left, right, pindex)


pval = A(pindex); 
A = swap(A, pindex, right); 

storeindex = left; 

for i =left:right-1
    if A(i)>pval
        A =swap(A, storeindex, i);
        storeindex = storeindex+1;
    end
     
end
A = swap(A, right, storeindex); 

end
function pivot = medofmed(A, left, right)
%for elements of 5 or less just get the median
if right - left <5
    pivot =floor(median(A(left:right)));
else
    for i = left:5:right
        subright = i+4; 
        if subright > right
            subright = right;
        end
        ind = i:subright; 
        median5 = floor(median(A(ind)));
%         median5 = ind(A(ind)==median5);  
        A = swap(A, median5, left+floor((i-left)/5));
        
    end
    pivot = select(A, left, left+floor((right-left)/5), (right-left)/10+1);
    
end


end
function A = swap(A, k1, k2)
%swap element of k1 with k2
if k1~=k2
    temp = A(k1); 
    A(k1) = A(k2); 
    A(k2) = temp; 
end
end