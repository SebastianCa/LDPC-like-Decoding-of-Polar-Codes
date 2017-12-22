%Naive dense polar H matrix according to Lemma 1 here: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5592698
%27.11.2017, Sebastian Cammerer, Moustafa Ebada, Ahmed Elkelesh, Stephan
%ten Brink
%{cammerer,ebada,elkelesh,tenbrink}@inue.uni-stuttgart.de

function [H, G_true] = createDensePolarH(N,A)
A=bitrevorder(A);   %legacy support for the rest of the framework
G=bitrevorder(computeG_method2(N)); 
G_true=G(A,:);
H=G(:,~A)';

if sum(sum(mod(H*G_true',2)))   %should be zero
    display('warning - something seems to be broken H*G!=0')
end

H=sparse(H);

end

function [G] = computeG_method2(N)
% This function calculates the Generator matrix for the polar code in the
% way mentioned in paper: Systematic polar codes, Erdal Arikan,  equation (6)
% A recursive function
% Inputs N Current code length
% Output G Generator matrix
if N == 1
    G =1;
elseif N==2 % 1st Stopping Condition
    G=[1 0; 1 1];
else
    G = ([ computeG_method2(N/2) zeros(N/2,N/2); computeG_method2(N/2) computeG_method2(N/2)  ]);
end
end
