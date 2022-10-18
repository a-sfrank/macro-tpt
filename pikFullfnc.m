function [pik_full]=pikFullfnc(k,chi_crisp,Q,target,n)
%% Description: 
% Function calculates the density peaks from the clustering result from the
% PCCA+ analysis.
%% Implemented by Susanna Röblitz (susanna.roblitz@uib.no)
        cluster=find(chi_crisp(:,k)==1);
        Bk=Q(cluster,cluster);
        [pik,~]=eigs(Bk,1,target);
        pik=pik/sum(pik);
        pik_full=zeros(n,1);
        pik_full(cluster)=pik;
end