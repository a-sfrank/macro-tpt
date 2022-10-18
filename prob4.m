function []=prob4(pi,Q,setA,setB,setC,Pc,qp,chi_crisp,fAB,setB1,setB2,tau)
%% Tristable case (Parameter Case 5, Subcases 3 to 6)
% Code implements Transition path theory (TPT) by
% (i) Calcuating stationary probability of a cluster of states (piTilde) (based on article (1));
% (ii) Transition probabilities between metastable sets (TIJTilde) (based on article (1));
% (iii) Transition flow betwenn metastable states to determine probabilistic transition paths (based on article (2));
%--------------------------------------
%% based on the articles:
% (1) Chu, Brian K., J. Tse Margaret, Royce R. Sato, and Elizabeth L. Read.
% "Markov State Models of gene regulatory networks." BMC systems biology 11, no. 1 (2017): 1-17.
% (2) Noé, Frank, Christof Schütte, Eric Vanden-Eijnden, Lothar Reich, and Thomas R. Weikl. 
% "Constructing the equilibrium ensemble of folding pathways from short off-equilibrium simulations."
%  Proceedings of the National Academy of Sciences 106, no. 45 (2009): 19011-19016.
%% Implemented by Anna-Simone Frank (anna-simone.frank@uib.no)
%--------------------------------------
% Code calcuates stationary probability of a cluster of states piTilde;
% Transition probabilities between metastable states TIJTilde
%--------------------------------------
%% (i) Calcuating stationary probability of a cluster of states
disp('Stationary probability of states:')
%Total stationary probability of a metastable states A and B:
piTildeA=sum(pi(setA)); X = sprintf('piA = %d',piTildeA); disp(X)
piTildeB=sum(pi(setB)); X = sprintf('piB = %d',piTildeB); disp(X)
piTildeB1=sum(pi(setB1)); X = sprintf('piB1 = %d',piTildeB1); disp(X)
piTildeB2=sum(pi(setB2)); X = sprintf('piB2 = %d',piTildeB2); disp(X)
piTildeC=sum(pi(setC)); X = sprintf('piC = %d',piTildeC); disp(X)

%------
% Coarse grained stationary denstiy --- Note: not used for article results.
% piTilde=pi.*chi_crisp;
%------
%-------------------------------------------------------------------------
% %% (ii) Transition probabilities between metastable sets 
% Note: not used for the results in the article.

% % Full transition probabilitiy:
% tau=100;
% Pfull=expm(tau.*Q);
% 
% % Transition prob within sets:
% disp('Transition probabilities within states:')
% TAATilde=sum(sum(Pfull(setA,setA).*pi(setA)))./piTildeA;  X = sprintf('TAA = %d',TAATilde); disp(X)
% TBBTilde=sum(sum(Pfull(setB,setB).*pi(setB)))./piTildeB; X = sprintf('TBB = %d',TBBTilde); disp(X)
% TB1B1Tilde=sum(sum(Pfull(setB1,setB1).*pi(setB1)))./piTildeB1; X = sprintf('TB1B1 = %d',TB1B1Tilde); disp(X)
% TB2B2Tilde=sum(sum(Pfull(setB2,setB2).*pi(setB2)))./piTildeB2; X = sprintf('TB2B2 = %d',TB2B2Tilde); disp(X)
% TCCTilde=sum(sum(Pfull(setC,setC).*pi(setC)))./piTildeC; X = sprintf('TCC = %d',TCCTilde); disp(X)
% 
% %Transition prob between sets:
% disp('Transition probabilities between states:')
% TABTilde=sum(sum(Pfull(setA,setB).*pi(setA)))./piTildeA; X = sprintf('TAB = %d',TABTilde); disp(X)
% TBATilde=sum(sum(Pfull(setB,setA).*pi(setB)))./piTildeB; X = sprintf('TBA = %d',TBATilde); disp(X)
% TB1ATilde=sum(sum(Pfull(setB1,setA).*pi(setB1)))./piTildeB1; X = sprintf('TB1A = %d',TB1ATilde); disp(X)
% TB2ATilde=sum(sum(Pfull(setB2,setA).*pi(setB2)))./piTildeB2; X = sprintf('TB2A = %d',TB2ATilde); disp(X)
% TB1B2Tilde=sum(sum(Pfull(setB1,setB2).*pi(setB1)))./piTildeB1; X = sprintf('TB1B2 = %d',TB1B2Tilde); disp(X)
% TB2B1Tilde=sum(sum(Pfull(setB2,setB1).*pi(setB2)))./piTildeB2; X = sprintf('TB2B1 = %d',TB2B1Tilde); disp(X)
% 
% TAB1Tilde=sum(sum(Pfull(setA,setB1).*pi(setA)))./piTildeA; X = sprintf('TAB1 = %d',TAB1Tilde); disp(X)
% TAB2Tilde=sum(sum(Pfull(setA,setB2).*pi(setA)))./piTildeA; X = sprintf('TAB2 = %d',TAB2Tilde); disp(X)
% TACTilde=sum(sum(Pfull(setA,setC).*pi(setA)))./piTildeA; X = sprintf('TAC = %d',TACTilde); disp(X)
% TCBTilde=sum(sum(Pfull(setC,setB).*pi(setC)))./piTildeC; X = sprintf('TCB = %d',TCBTilde); disp(X)
% TBCTilde=sum(sum(Pfull(setB,setC).*pi(setB)))./piTildeB; X = sprintf('TBC = %d',TBCTilde); disp(X)
% TCB1Tilde=sum(sum(Pfull(setC,setB1).*pi(setC)))./piTildeC; X = sprintf('TCB1 = %d',TCB1Tilde); disp(X)
% TCATilde=sum(sum(Pfull(setC,setA).*pi(setC)))./piTildeC; X = sprintf('TCA = %d',TCATilde); disp(X)
% TCB2Tilde=sum(sum(Pfull(setC,setB2).*pi(setC)))./piTildeC; X = sprintf('TCB2 = %d',TCB2Tilde); disp(X)
% 
% TB1CTilde=sum(sum(Pfull(setB1,setC).*pi(setB1)))./piTildeB1; X = sprintf('TB1C = %d',TB1CTilde); disp(X)
% TACTilde=sum(sum(Pfull(setA,setC).*pi(setA)))./piTildeA; X = sprintf('TAC = %d',TACTilde); disp(X)
% TB2CTilde=sum(sum(Pfull(setB2,setC).*pi(setB2)))./piTildeB2; X = sprintf('TB2C = %d',TB2CTilde); disp(X)
%-------------------------------------------------------------
%% (iii) Transition flow betwenn metastable states to determine probabilistic transition paths
% Total amount of flux out of A and into B:
disp('Total transition flux from A and into B:')
setBC=union(setC,setB);
setAC=union(setA,setC);
Ftot2a=sum(fAB(setA,setBC),'all'); X = sprintf('Flux out of A: %d',Ftot2a); disp(X)
Ftot2b=sum(fAB(setAC,setB),'all'); X = sprintf('Flux into B: %d',Ftot2b); disp(X)


% Flow decomposition for coarse grained flux:
disp('Flow decompositions for coarse grained flux:')
% (From this decomposed flux, one can caluclate the relative
% probabilities.)
% Redefine C, s.t., C=C\B2. For the flow we decompose trapped state from
% set C:
setC=setdiff(setC,setB2); 
%disp('Flow decompositions for coarse grained flux:')
Ftot3a1=sum(fAB(setA,setB1),'all'); X = sprintf('Flux A->B1: %d',Ftot3a1); disp(X)
Ftot3b1=sum(fAB(setA,setB2),'all'); X = sprintf('Flux A->B2: %d',Ftot3b1); disp(X) 
Ftot3c1=sum(fAB(setA,setC),'all'); X = sprintf('Flux A->C: %d',Ftot3c1); disp(X) 
Ftot3d1=sum(fAB(setA,setA),'all'); X = sprintf('Flux A->A: %d',Ftot3d1); disp(X) 

Ftot4a1=sum(fAB(setB1,setA),'all'); X = sprintf('Flux B1->A: %d',Ftot4a1); disp(X) 
Ftot4b1=sum(fAB(setB1,setB2),'all'); X = sprintf('Flux B1->B2: %d',Ftot4b1); disp(X)  
Ftot4c1=sum(fAB(setB1,setC),'all'); X = sprintf('Flux B1->C: %d',Ftot4c1); disp(X)  

Ftot5a1=sum(fAB(setB2,setA),'all'); X = sprintf('Flux B2->A: %d',Ftot5a1); disp(X)  
Ftot5b1=sum(fAB(setB2,setB1),'all'); X = sprintf('Flux B2->B1: %d',Ftot5b1); disp(X)  
Ftot5c1=sum(fAB(setB2,setC),'all'); X = sprintf('Flux B2->C: %d',Ftot5c1); disp(X) 
Ftot5d1=sum(fAB(setB2,setB2),'all'); X = sprintf('Flux B2->B2: %d',Ftot5d1); disp(X) 

Ftot6a1=sum(fAB(setC,setA),'all'); X = sprintf('Flux C->A: %d',Ftot6a1); disp(X)  
Ftot6b1=sum(fAB(setC,setB1),'all'); X = sprintf('Flux C->B1: %d',Ftot6b1); disp(X)  
Ftot6c1=sum(fAB(setC,setB2),'all'); X = sprintf('Flux C->B2: %d',Ftot6c1); disp(X)  