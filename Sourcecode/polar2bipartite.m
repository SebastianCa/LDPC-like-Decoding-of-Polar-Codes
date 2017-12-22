%BP Decoding of Polar Codes using Matlab's conventional BP decoder
%see paper "Sparse Graphs for Belief Propagation Decoding of Polar Codes"
%27.11.2017, Sebastian Cammerer, Moustafa Ebada, Ahmed Elkelesh, Stephan
%ten Brink
%{cammerer,ebada,elkelesh,tenbrink}@inue.uni-stuttgart.de


function [H]=polar2bipartite(N,A,permute)

if nargin<3
    permute=0;
end
n=log2(N);  %number of stages

nbCN=N;
nbVN=N*(n+1);

%initalization
nzmax=3*nbVN;       %number of nonzero elements in H

H = spalloc(nbCN,nbVN,nzmax);   %initlaize sparse H matrix
CNid=1;

%and let´s go
I = 1:n;

%permute I here, see WCNC 2018 "BP decoding of polar codes on permuted factor graphs"
if permute
    permInd=randperm(n-2);
    t=2:(n-1);
    I(2:(n-1))=t(permInd);
end
for id=1:length(I) % log2(N) stages (n stages)
    i = I(id);
    temp=-1*ones(1,N); % In order not to visit a bit twice
    divisor=2^i;
    separation=N/divisor; % Separation between the bits to be added together
    for j=1:N % To pass through the N bits at every stage
        if temp(j)==-1 % In order not to visit a bit twice
            temp(j)=0; % Bit at index j visited (output at index j for the current stage calculated already)
            temp(j+separation)=0; % Bit at index j+separation visited (output at index j+separation for the current stage calculated already)
            
            vn1=j+(i-1)*N;
            vn2=j+separation+(i-1)*N;
            
            %forward connection
            H(CNid,vn1)=1;
            H(CNid,vn2)=1;
            %backward connection from next layer
            H(CNid,vn1+N)=1;
            %next CN
            CNid=CNid+1;
            H(CNid,vn2)=1;
            %backward connection from next layer
            H(CNid,vn2+N)=1;
            CNid=CNid+1;
            
        end
    end
end

%kick out VNs with frozen indices
A_nonfrozen=ones(size(H,2),1);
A_nonfrozen(1:N)=A;
A_nonfrozen=logical(A_nonfrozen);
H=H(:,A_nonfrozen);

end

