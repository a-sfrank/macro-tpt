function []=prob2(pi,Q,setA,setB,setC,Pc,qp,chi_crisp,fAB,tau)
%% Bistable cases (Parameter Cases 1 to4)
% Code implements Transition path theory (TPT) 
% (i) Calcuating stationary probability of a cluster of states (piTilde) (based on article (1));
% (ii) Transition probabilities between metastable sets (TIJTilde) (based on article (1));
% (iii) Transition flow betwenn metastable states to determine probabilistic transition paths (based on article (2));
%--------------------------------------
%% based on articles:
% (1) Chu, Brian K., J. Tse Margaret, Royce R. Sato, and Elizabeth L. Read.
% "Markov State Models of gene regulatory networks." BMC systems biology 11, no. 1 (2017): 1-17.
% (2) Noé, Frank, Christof Schütte, Eric Vanden-Eijnden, Lothar Reich, and Thomas R. Weikl. 
% "Constructing the equilibrium ensemble of folding pathways from short off-equilibrium simulations."
%  Proceedings of the National Academy of Sciences 106, no. 45 (2009): 19011-19016.
%% Implemented by Anna-Simone Frank (anna-simone.frank@uib.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------
%% (i) Total stationary probability of a metastable states A and B:
disp('Stationary probability of states:')
piTildeA=sum(pi(setA)); X = sprintf('piA = %d',piTildeA); disp(X)
piTildeB=sum(pi(setB)); X = sprintf('piB = %d',piTildeB); disp(X)
piTildeC=sum(pi(setC)); X = sprintf('piC = %d',piTildeC); disp(X)

%-------------------------------------------------------------------------
%% (ii) Transition probabilities between metastable sets: 
% Note: Commented out because not used in the article;

% Full transition probabilitiy
% Pfull=expm(tau.*Q);
% 
% % Transition prob within metastable sets:
% disp('Transition probabilities within states:')
% TAATilde=sum(sum(Pfull(setA,setA).*pi(setA)))./piTildeA; X = sprintf('TAA = %d',TAATilde); disp(X)
% TBBTilde=sum(sum(Pfull(setB,setB).*pi(setB)))./piTildeB; X = sprintf('TBB = %d',TBBTilde); disp(X)
% TCCTilde=sum(sum(Pfull(setC,setC).*pi(setC)))./piTildeC; X = sprintf('TCC = %d',TCCTilde); disp(X)
% 
% %Transition prob between sets:
% disp('Transition probabilities between states:')
% TABTilde=sum(sum(Pfull(setA,setB).*pi(setA)))./piTildeA; X = sprintf('TAB = %d',TABTilde); disp(X)
% TBATilde=sum(sum(Pfull(setB,setA).*pi(setB)))./piTildeB; X = sprintf('TBA = %d',TBATilde); disp(X)
% TACTilde=sum(sum(Pfull(setA,setC).*pi(setA)))./piTildeA; X = sprintf('TAC = %d',TACTilde); disp(X)
% TCBTilde=sum(sum(Pfull(setC,setB).*pi(setC)))./piTildeC; X = sprintf('TCB = %d',TCBTilde); disp(X)
% TBCTilde=sum(sum(Pfull(setB,setC).*pi(setB)))./piTildeB; X = sprintf('TBC = %d',TBCTilde); disp(X)
% TCATilde=sum(sum(Pfull(setC,setA).*pi(setC)))./piTildeC; X = sprintf('TCA = %d',TCATilde); disp(X)

%--------------------------------------------------------------------------
%% (iii) Transition flow betwenn metastable states to determine probabilistic transition paths 
% Total amount of flux out of A and into B:
disp('Total transition flux from A and into B:')
setBC=union(setC,setB);
setAC=union(setA,setC);
Ftot2a=sum(fAB(setA,setBC),'all'); X = sprintf('Flux out of A: %d',Ftot2a); disp(X)
Ftot2b=sum(fAB(setAC,setB),'all');  X = sprintf('Flux into B: %d',Ftot2b); disp(X)

% Flow decomposition for coarse grained flux:
disp('Flow decompositions for coarse grained flux:')
% (From this decomposed flux, one can caluclate the relative
% probabilities )
Ftot3a1=sum(fAB(setA,setB),'all'); X = sprintf('Flux A->B: %d',Ftot3a1); disp(X)
Ftot3b1=sum(fAB(setA,setC),'all'); X = sprintf('Flux A->C: %d',Ftot3b1); disp(X)

Ftot4a1=sum(fAB(setB,setA),'all'); X = sprintf('Flux B->A: %d',Ftot4a1); disp(X)
Ftot4b1=sum(fAB(setB,setC),'all'); X = sprintf('Flux B->C: %d',Ftot4b1); disp(X)

Ftot6a1=sum(fAB(setC,setA),'all'); X = sprintf('Flux C->A: %d',Ftot6a1); disp(X)
Ftot6b1=sum(fAB(setC,setB),'all'); X = sprintf('Flux C->B: %d',Ftot6b1); disp(X)