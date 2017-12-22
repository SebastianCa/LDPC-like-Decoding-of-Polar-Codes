%create the Polar code curves for the paper "Sparse Graphs for Belief Propagation Decoding of Polar Codes"
%27.11.2017, Sebastian Cammerer, Moustafa Ebada, Ahmed Elkelesh, Stephan
%ten Brink
%{cammerer,ebada,elkelesh,tenbrink}@inue.uni-stuttgart.de


function createPaperPlots()

%general parameters
simparam.R=0.5;
simparam.desSNR=0.6;
simparam.nbSNR=15;
simparam.useSysEncoding=1;   %use systematic encoding
simparam.bipartite=1;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)

%specific parameters
simparam.N=256;
simparam.iterMax=200;
simparam.nbCWs=1e5;
simparam.SNRb_start=0;
simparam.SNRb_stop=6;
simparam.bipartite=1;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)
simparam.filename='N256-R05-SNR06';

startSim(simparam);

simparam.SNRb_start=0;
simparam.SNRb_stop=10;
simparam.bipartite=0;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)
simparam.filename='N256-R05-SNR06-Dense';

startSim(simparam);


%specific parameters
simparam.N=2048;
simparam.iterMax=200;
simparam.nbCWs=1e5;
simparam.SNRb_start=0;
simparam.SNRb_stop=5;
simparam.bipartite=1;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)
simparam.filename='N2048-R05-SNR06';

startSim(simparam);

simparam.SNRb_start=0;
simparam.SNRb_stop=10;
simparam.bipartite=0;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)
simparam.filename='N2048-R05-SNR06-Dense';

startSim(simparam);


%specific parameters
simparam.N=32768;
simparam.iterMax=200;
simparam.nbCWs=1e5;
simparam.SNRb_start=0;
simparam.SNRb_stop=4;
simparam.bipartite=1;        %use LDPC-like code construction as presented in the paper (naive dense matrix is used instead)
simparam.filename='N32768-R05-SNR06';

startSim(simparam);

end


function startSim(simparam)
simres=runSim(simparam);  %run simulation
save(['results/' simparam.filename],'simparam','simres');    %save results
export_pgfplots( simres.SNRb, simres.BER, ['results/' simparam.filename '-BER']);   %save BER curves for pgfplots import
export_pgfplots( simres.SNRb, simres.BER, ['results/' simparam.filename '-BLER']);   %save BLER curves for pgfplots import
end