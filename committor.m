function [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q,pi,setA,setB,Case,filepath)
% Code implements commitor functions from article: (1) Metzner, Philipp, Christof Schütte, and Eric Vanden-Eijnden. 
% "Transition path theory for Markov jump processes." Multiscale Modeling & Simulation 7.3 (2009): 1192-1219.
%% Implemented by Susanna Röblitz (susanna.roblitz@uib.no)
% Notes:
%Q: rate matrix with row sum 0
%set A: indices of states belonging to set A
%set B: indices of states belonging to set B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n=size(Q,1);
    N=sqrt(n);
    
    setC=setdiff([1:n],setA);
    setC=setdiff(setC,setB);
    setC=setC';
       
    qpC=Q(setC,setC)\(-sum(Q(setC,setB),2)); %eq. (2.13.) from (1)
    qp=zeros(n,1);
    qp(setC)=qpC;
    qp(setA)=0;
    qp(setB)=1;
    figure(21)
    h=surf(reshape(qp,N,N));
    set(h,'LineStyle','none')
    title('forward committor') 
    fname= sprintf('Fig_21_Case_%d.png',Case);
    saveas(h, fullfile(filepath,fname));

    Qtilde=diag(pi)*Q*diag(1./pi);  %eq. (2.4) from (1)
    Qtilde=Qtilde';

    qmC=Qtilde(setC,setC)\(-sum(Q(setC,setA),2));      %eq. (2.14.) from (1)
    qm=zeros(n,1);
    qm(setC)=qmC;
    qm(setA)=1;
    qm(setB)=0;
    figure(22)
    h=surf(reshape(qm,N,N));
    set(h,'LineStyle','none')
    title('backward committor')
    fname= sprintf('Fig_22_Case_%d.png',Case);
    saveas(h, fullfile(filepath,fname));
    
    %filename='ComittorVars.mat';
    %save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n')

end