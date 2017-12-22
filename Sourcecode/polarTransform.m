function [x] = polarTransform(u, A)
% non-systematic polar encoding
N=length(u);

for i=1:log2(N) % log2(N) stages (n stages)
    temp=-1*ones(1,N); % In order not to visit a bit twice
    divisor=2^i;
    separation=N/divisor; % Separation between the bits to be added together
    for j=1:N % To pass through the N bits at every stage
        if temp(j)==-1 % In order not to visit a bit twice
            u(j)=mod(u(j)+u(j+separation),2); % First output of the butterfly, the second output will not change (equal to the second input)
            temp(j)=0; % Bit at index j visited (output at index j for the current stage calculated already)
            temp(j+separation)=0; % Bit at index j+separation visited (output at index j+separation for the current stage calculated already)
        end
    end
end
x=u;

end