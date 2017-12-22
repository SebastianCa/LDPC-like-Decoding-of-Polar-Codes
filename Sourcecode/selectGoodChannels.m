function A = selectGoodChannels(N, k, SNR)
% INPUT : N, K, and SNR design-SNR  = (REb/N0 in dB)
% OUTPUT: A subset of {0, 1, . . . , N - 1} with |A==1| = k logical [1 X N]

Ac = false(1,N); % initially choose No channels
S = 10^(SNR/10); 
n = log2(N);
z0(1) = exp(-S);
for j = 1:n
    u = 2^j;
    for t = 1:u/2
        T = z0(t);
        z0(t) = 2*T - T^2   ; % Upper Bad channel
        z0(t + u/2) = T^2   ; % Lower Good Channel
    end
end

Ac(FindMaxElem(z0, N-k)) = 1; % Find indices of the greatest N - K elements
A = bitrevorder(~Ac);
end
function I = FindMaxElem(v, k)
% INPUT v: vector of dimension |v| W 1 
%       k: integer 
% OUTPUT I: k W 1 integer vector containing k indices in {0, 1, . . . , |v| - 1}
[~, idx] = sort(v, 'descend'); % obtain in idx, the |v| indices of vector v when sorted in descending order
I = idx(1:k); % B Store the first k indices
end