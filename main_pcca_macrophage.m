%% Main run file to perform analysis presented in the article "Macrophage phenotype transitions in a stochastic 
%% gene-regulatory netwrok model" by Anna S. Frank, Kamila Larripa, Hwayeon Ryu and Susanna Röblitz.
% Code loads the stochastic macrophage model, performs the PCCA+ and TPT analysis 
%% Implemented by Susanna Röblitz (susanna.roblitz@uib.no), Anna S. Frank (anna-simone.frank@uib.no)
close all; clear all;
tic
% load pre-specified information:
Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary('.\Results_Main_PCCA\PCCA_results_File.txt')
diary(filePCCA)
diary on
%% evtl. Time and date stamp or if exist then remove file ????

%generate the CME matrix (column sums are zero)
 [Q,Case]=assembleQ(N);  % assembleQ file calls parameter cases
 assert(max(abs(sum(Q,1)))<1e-12,'Columns sums of Q are not zero!')
 n=size(Q,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Coarse-graining and clustering via PCCA+
%density for normalization of eigenvectors
 ovec=ones(n,1)/n; 

% Using P:
   P=expm(tau*Q);
   [n,m]=size(P);

% Compute Sub-space
[EVS,la,kmax]=compute_subspace(P',kmax,target,opts); 
Y=sprintf('Real part of eigenvalues: la'); 
disp(Y)
la

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization step by PCCA+: Compute A such that chi=EVS*A 
[chi,A,optval,EVS]=pcca(EVS,ovec,kmin,kmax,flag,solver,Case,filepath);

%number of clusters:
nc=size(chi,2);

disp (' ')
disp ('Coarse-grainded transition matrix:') 

% Calculating Qc:
Qc=inv(chi'*diag(ovec)*chi)*chi'*diag(ovec)*Q'*chi  


%Using P: Calculate coarse-grained transition matrix Pc
Pc=inv(chi'*diag(ovec)*chi)*chi'*diag(ovec)*P'*chi  


disp('===================================================================')
diary off

% Plot the membership functions in 2D
for k=1:nc
    chi_2D=reshape(chi(:,k),N+1,N+1);
    figure(k+1)
    h=contour(chi_2D);
    colorbar
    title(sprintf('Membership functions, chi=%d',k))
    fname1= sprintf('Fig_%d_Case_%d.png', k+1, Case);
    saveas(gcf, fullfile(filepath,fname1));
     fname2= sprintf('Fig_%d_Case_%d.eps', k+1, Case);
    saveas(gcf, fullfile(filepath,fname2));
     fname3= sprintf('Fig_%d_Case_%d.fig', k+1, Case);
    saveas(gcf, fullfile(filepath,fname3));
end


% Make the fuzzy clustering crisp
chi_crisp=zeros(n,nc);
for i=1:n
    [val,idx]=max(chi(i,:));
    chi_crisp(i,idx)=1;
end

pik_2D_sum=zeros(N+1,N+1);
for k=1:nc
        [pik_full]=pikFullfnc(k,chi_crisp,Q,target,n);
        findSet_2(k,pik_full,filepath,Case,N)
        %------------------------------------------------------------------
        pik_2D=reshape(pik_full,N+1,N+1);
        pik_2D_sum=pik_2D_sum+pik_2D;
        figure(10+k)
        h=surf(pik_2D);
        set(h,'LineStyle','none');
        colorbar
        title(sprintf('Partial stationary densities of clusters, k=%d', k))
        fname= sprintf('Fig_%d_Case_%d.png', 10+k, Case);
        saveas(h, fullfile(filepath,fname));
        fname2= sprintf('Fig_%d_Case_%d.eps', 10+k, Case);
        saveas(h, fullfile(filepath,fname2));
        fname3= sprintf('Fig_%d_Case_%d.fig', 10+k, Case);
        saveas(h, fullfile(filepath,fname3));
end
figure(10+nc+1)
h=surf(pik_2D_sum);
set(h,'LineStyle','none');
colorbar
set(gca,'ColorScale','log')
title(sprintf('Partial stationary densities'));
fname= sprintf('Fig_%d_Case_%d.png', 10+nc+1, Case);
saveas(gcf, fullfile(filepath,fname));
fname2= sprintf('Fig_%d_Case_%d.eps', 10+nc+1, Case);
saveas(gcf, fullfile(filepath,fname2));
fname3= sprintf('Fig_%d_Case_%d.fig', 10+nc+1, Case);
saveas(gcf, fullfile(filepath,fname3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                TPT analysis
[pi,la]=eigs(Q,1,1+eps); % Eigenvalues from rate matrix Q are calculated
[val,idx]=find(max(abs(pi)));
pi=pi*sign(pi(idx));
pi=max(pi,eps);
pi=pi/sum(pi);

TPTCases(Case,Q,pi,filepath,Pc,chi_crisp,tau)



    




