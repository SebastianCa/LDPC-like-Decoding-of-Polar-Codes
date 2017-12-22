%BP Decoding of Polar Codes using Matlabs conventional BP decoder
%see paper "Sparse Graphs for Belief Propagation Decoding of Polar Codes"
%27.11.2017, Sebastian Cammerer, Moustafa Ebada, Ahmed Elkelesh, Stephan
%ten Brink
%{cammerer,ebada,elkelesh,tenbrink}@inue.uni-stuttgart.de

%   master pruning function (Algorithm 1), removes all unneeded VNs and CNs

function [ H ] = pruneGraph( H,N )

size_old=Inf;
while size_old>sum(size(H))    %repeat until size(H) is does not change
    size_old=sum(size(H));
    
    H=pruneDegree1CN( H );  %prune all CN of degree and delete connected VNs, as they are forced to zero anyway
    
    H=condenseDegree1VNch( H,N );
    
    H=pruneDegree1VNH(H,N);
    
    H=condenseDegree2VNH(H,N);      %remove all hidden degree 2 VNs as they simply forward messages
    
    H=condenseDegree2CN(H,N);
    
    H=deleteEmptyChecks(H,N);
    
end

end


function [ H ] = pruneDegree1CN( H )
%identify CN of degree 1 and kick out
%also kick out the connected VN (as they are forced to 0 anyway)

prune=1;
while prune
    CN1_ind=find(sum(H,2)==1);
    Htemp=H(CN1_ind,:);
    Avn=sum(Htemp,1)>0;
    H(CN1_ind,:)=[];        %kickout CN        
    H(:,Avn)=[];
    if length(CN1_ind)==0   %prune until no CN-1 remains
        prune=0;
    end
end

end

function [ H ] = condenseDegree1VNch( H,N )
%   if a channel VN is of degree 1 and connected to a CN of degree 2, then it only propagates Lch to the next VN
%idea: kick out this VN and assign Lch to connected VN (i.e., swap columns)

for i=0:N-1
   c=H(:,end-i);
   if sum(c)==1 %if VN is of degree 1
       ind=find(c);
       r=H(ind,:);
       if sum(r)==2
          ind2=find(r);
          ind2=ind2(1);
          H(:,end-i)=H(:,ind2);
          H(:,ind2)=c;
       end
   end           
end
%call pruneDegree1VNH later to kick out the swaped nodes
end


function [ H ] = pruneDegree1VNH( H,N )
%idea: each CN which is connected to at least one hiddenVN of degree 1 only
%contributes 0 (or close to 0) to all VNs => remove the CN

keepPruning=1;
while keepPruning
    size_prev=sum(size(H));
    nbVN=size(H,2);
    nbCN=size(H,1);
    %idea: each CN which is connected to at least one hiddenVN of degree 1 only
    %contributes 0 (or close to 0) to all VNs => remove the CN
    A=ones(nbCN,1);
    for i=1:nbCN
        r=H(i,:);
        ind=find(r(1:end-N));    %only hidden positions are considered
        for j=1:length(ind)
            c=H(:,ind(j));
            if sum(c)==1  %remove cn
                A(i)=0;
            end
        end        
    end
    %slice H
    H=H(logical(A),:);
    %remove all-zero columns
    H(:,(sum(H,1)==0))=[];
    if size_prev==sum(size(H));
        keepPruning=0;
    end
end


end


function [ H ] = condenseDegree2VNH( H,N )
%   idea: degree 2 hidden VNs are just forward nodes (as they never have
%   own internal information). In case they are connected to a CN of degree
%   2, there can be a short cut

nbVN=size(H,2);
i=1;
while i<nbVN-N+1
    c=H(:,i);
    if sum(c)==2
        ind=find(c);
        ind1=ind(1);
        ind2=ind(2);
        H(ind1,:)=mod(H(ind1,:)+H(ind2,:),2);
        H(ind2,:)=[];
        H(:,i)=[];
    end
    nbVN=size(H,2);
    i=i+1;
end

end

 function [ H ] = condenseDegree2CN( H,N )
%   idea a degree 2 CN simply forwards information, therefore we can build
%   a "super hidden" node which represents both...do not do this for
%   channel nodes

nbCN=size(H,1);
i=1;
while i<nbCN
    r=H(i,:);
    if sum(r(1:end-N))==2 && sum(r(end-N+1:end))==0
        ind=find(r);
        ind1=ind(1);
        ind2=ind(2);
        H(:,ind1)=mod(H(:,ind1)+H(:,ind2),2);
        H(:,ind2)=[];
    end
    nbCN=size(H,1);
    i=i+1;
end

 end

 function [ H ] = deleteEmptyChecks(H,N)
% delete empty CN and VNs as they do nothing (except for channel nodes in
% order to match the dimension requirements
nbVN=size(H,2);
nbCN=size(H,1);

A=ones(nbCN,1);
for i=1:nbCN
   r=H(i,:);
   if sum(r)==0
       A(i)=0;
   end
end
H=H(logical(A),:);

A=ones(nbVN,1);
for i=1:nbVN-N;
   c=H(:,i);
   if sum(c)==0
       A(i)=0;
   end
end
H=H(:,logical(A));

end
