%BP Decoding of Polar Codes using Matlabs conventional BP decoder
%see paper "Sparse Graphs for Belief Propagation Decoding of Polar Codes"
%27.11.2017, Sebastian Cammerer, Moustafa Ebada, Ahmed Elkelesh, Stephan
%ten Brink
%{cammerer,ebada,elkelesh,tenbrink}@inue.uni-stuttgart.de

function [simres]=runSim(simparam)

%Simulation Parameters
if nargin==0
    simparam.N=256;
    simparam.R=0.5;
    simparam.desSNR=0.6;
    simparam.iterMax=200;
    simparam.nbCWs=1e3;
    simparam.nbSNR=4;
    simparam.SNRb_start=1;
    simparam.SNRb_stop=4;
    simparam.useSysEncoding=1;   %use systematic encoding
    simparam.bipartite=1;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)
end


%Optional features

useAllZeroCW=0;     %use the all zero CW for the simulations
createPlots=0;      %plot results

%Initialize variables and
simres.SNRb=linspace(simparam.SNRb_start,simparam.SNRb_stop,simparam.nbSNR);
simparam.k=round(simparam.N*simparam.R);

simres.nbBits=zeros(simparam.nbSNR,1);
simres.nbErr=zeros(simparam.nbSNR,1);
simres.nbBlocks=zeros(simparam.nbSNR,1);
simres.nbBlockErr=zeros(simparam.nbSNR,1);

%create information set
A=selectGoodChannels(simparam.N, simparam.k, simparam.desSNR);
A=bitrevorder(A);

%%%%%%%%%start pruning%%%%%%%%%%%%%%5
display('Create and prune H')
%create H matrix

if simparam.bipartite==1
    [H]=polar2bipartite(simparam.N,A,1); %already kicks out frozen indices
    size(H)    
    %prune all useless nodes
    H=pruneGraph(H,simparam.N);
    size(H)
    simres.A=A;
    simres.Hpruned=H;    
else
    H=createDensePolarH(simparam.N,A);
    simres.Hdense=H;
    size(H)    
    spy(H)    
end

if createPlots==1
    figure(1);clf;hold;
    spy(H)  %show H matrix
end

%%%%%%%%%%%%START BER SIMULATIONS%%%%%%%%%%%%%%%%%
%initialize decoder
display('Starting BER simulation')

decoder = comm.LDPCDecoder('ParityCheckMatrix',H,'OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount',simparam.iterMax,'NumIterationsOutputPort',1,'IterationTerminationCondition', 'Parity check satisfied');

for j=1:simparam.nbSNR
    display(['evaluating SNRb=' num2str(simres.SNRb(j)) 'dB, ' num2str(j) '/' num2str(simparam.nbSNR)]);
    
    sigma=sqrt(0.5*((10^((simres.SNRb(j)+10*log10(simparam.R))/10))^-1));
    %initialize temp variables
    nbErr_temp=zeros(simparam.nbCWs,1);
    nbBit_temp=zeros(simparam.nbCWs,1);
    nbBlockErr_temp=zeros(simparam.nbCWs,1);
    nbBlocks_temp=zeros(simparam.nbCWs,1);
    numiter=zeros(simparam.nbCWs,1);
    
    parfor i=1:simparam.nbCWs
        
        if useAllZeroCW
            d=zeros(simparam.k,1); %create all-zero data
        else
            d=randi(2,simparam.k,1)-1; %create random data
        end
        u=zeros(simparam.N,1);
        u(A)=d;
        x=polarTransform(u, A);
        
        %reencode for systematic encoding, see Paper Arikan
        if simparam.useSysEncoding
            x(~A)=0;
            x=polarTransform(x, A);
        end
        
        x=-2*x+1;   %Bi-AWGN channel
        y=x+sigma*randn(simparam.N,1);
        Lch=2*y./sigma.^2;  %calc LLRs
        
        %extend Lch vector by 0 positions
        Lch_ext=zeros(size(H,2),1);
        Lch_ext((end-simparam.N+1):end)=Lch;     %assuming the last positions are channel positions
        
        [LLRxhat, iter] = step(decoder, Lch_ext);     %and decode
        numiter(i)=iter;    %save iteration count
        
        xhat=(-sign(LLRxhat))/2+0.5;
        
        c=xhat((end-simparam.N+1):end);
        if simparam.useSysEncoding==0
            c=polarTransform(c, A);%use x to predict c as this simplifies the notation in combination with pruning
        end
        c=c(A);
        
        %count errors
        nbErr_temp(i)=sum(c~=d);
        nbBit_temp(i)=length(c);
        %also consider BLER
        if sum(c~=d)~=0
            nbBlockErr_temp(i)=1;
        end
        nbBlocks_temp(i)=1;
    end
    %sum up temp arrays (due to parfor loop)
    simres.nbErr(j)=simres.nbErr(j)+sum(nbErr_temp);
    simres.nbBits(j)=simres.nbBits(j)+sum(nbBit_temp);
    simres.nbBlocks(j)=simres.nbBlocks(j)+sum(nbBlocks_temp);
    simres.nbBlockErr(j)=simres.nbBlockErr(j)+sum(nbBlockErr_temp);
    simres.iter(j)=mean(numiter);
end

%print results
simres.nbErr
simres.nbBits
simres.BER=simres.nbErr./simres.nbBits
simres.BLER=simres.nbBlockErr./simres.nbBlocks

if createPlots==1
    figure(2);clf;
    semilogy(simres.SNRb,BER);hold;
    semilogy(simres.SNRb,BLER,':');
    legend('BER','FER');title('BER for bipartite Polar decoding');xlabel('SNRb [dB]');ylabel('BER/FER');
end

end