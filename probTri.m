 function [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probTri(qp,qm,setC,pi,Q,setA,setB,n,N,Case,SubC,filepath)
% Code implements parts of articles: (1) Metzner, Philipp, Christof Schütte, and Eric Vanden-Eijnden. 
% "Transition path theory for Markov jump processes." Multiscale Modeling & Simulation 7.3 (2009): 1192-1219.
% (2) Helfmann, Luzie, Enric Ribera Borrell, Christof Schütte, and Péter Koltai.
% "Extending transition path theory: Periodically driven and finite-time dynamics." Journal of nonlinear science 30, no. 6 (2020): 3321-3366.
%% Implemented by Anna-Simone Frank (anna-simone.frank@uib.no) Susanna Röblitz (susanna.roblitz@uib.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin of TPT analysis
mR=zeros(n,1);
mR(setC)=pi(setC).*qp(setC).*qm(setC); % eq. (2.17) from (1)
mR(setA)=pi(setA).*qp(setA).*qm(setA);
mR(setB)=pi(setB).*qp(setB).*qm(setB);

figure(31)
h=surf(reshape(mR,N,N));
set(h,'LineStyle','none')
title('Probability distribution of reactive trajectories (mR)') 
fname= sprintf('Fig_31_Case_%d_SubC_%d.png',Case,SubC);
saveas(h, fullfile(filepath,fname));
fname2= sprintf('Fig_31_Case_%d_SubC_%d.eps',Case,SubC);
saveas(h, fullfile(filepath,fname2));
fname3= sprintf('Fig_31_Case_%d_SubC_%d.fig',Case,SubC);
saveas(h, fullfile(filepath,fname3));

ZAB=sum(mR,1) ;                           %eq.(2.18) from (1)
mAB=mR./ZAB;                              %eq.(2.20) from (1)
    
figure(32)
h=surf(reshape(mAB,N,N));
set(h,'LineStyle','none')
title('Normalized distribution of reactive trajectories (mAB)') 
fname= sprintf('Fig_32_Case_%d_SubC_%d.png',Case,SubC);
saveas(h, fullfile(filepath,fname));
fname2= sprintf('Fig_32_Case_%d_SubC_%d.eps',Case,SubC);
saveas(h, fullfile(filepath,fname2));
fname3= sprintf('Fig_32_Case_%d_SubC_%d.fig',Case,SubC);
saveas(h, fullfile(filepath,fname3));

% Reversibility
Vji=diag(pi)*Q*diag(1./pi);
T=isequal(Vji,Q');
if T==1
   disp("MC is reversible")
else
   disp("MC not reversible")
end
    
% Compute the whole state space
AB=union(setA,setB);
ABC=union(setC,AB);
n=length(ABC);
fAB=zeros(n,n);
size(Q);
[ii,jj]=find(Q);   %find indices of nonzero elements in sparse matrix Q

for k=1:length(ii)
    if ii(k)~=jj(k)
       fAB(ii(k),jj(k))=pi(ii(k)).*qm(ii(k)).*Q(ii(k),jj(k)).*qp(jj(k)) ; % eq. (2.24) from (1)
    end
end

% F effective
Feff=max(fAB(ABC,ABC)-fAB(ABC,ABC)',0);

fname=sprintf('case%d_SubCase%d_Feff.mat',Case,SubC);
save(fname,'Feff','setA','setB');


%% Plot figure 33 (optional) by uncommenting or run function separate:
% plot_Feff(Case,SubC)

%Conservation
C=sum(fAB(setC,ABC)-fAB(ABC,setC)',2);  %eq. (26)  from (1) %transpose
if max(abs(C))<10^(-6)
   disp('Conservation of discrete probability flux')
else
   disp('No conservation of flux')
end

%Reaction rate
kAB1=sum(fAB(setA,ABC),'all')  ;           %eq. (2.30) from (1)
kAB2=sum(fAB(ABC,setB),'all');
if (kAB1-kAB2)<10^(-6)
   disp('Transition rate as in eq.2.30')
else
   disp('error in calculation of transition rate')
end
 
%Expected length of a transition from A to B (from (2) Helfmann et al.)
tAB=ZAB./kAB1;

%save files (optional)
% filename='ProbVars.mat';
% save(filename,'mAB','fAB','Feff','kAB1','ZAB','tAB')
 

end
    
 
    
   
