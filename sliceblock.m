function Y = sliceblock(X,K,P,N)
%Slice repeated matrix blocks and rearrange them horizontally
%   X is a (K * N)* P matrix where the block matrices is K * P where are
%   repeated N times. Y would be a K * (P * N) matrix


Y = reshape(permute(reshape(X', P, K, N), [2 1 3]), K, N*P);

end

