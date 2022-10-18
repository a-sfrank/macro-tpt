function []=findSet_2(k,pik_full,filepath,Case,N) 
% Function code extracts the sets defined from coarse-graining from PCCA+
% clustering. These sets represent the phenotypes and are needed for the
% further TPT analysis. They are saved in .mat files and called in
% TPTCases.m file for the TPT analysis.
% A maximum of three (k=3) sets are considered, since the PCCA analysis
% resulted in a maximum of three sets.
%% Implemented by Anna S. Frank (anna-simone.frank@uib.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Input data
Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k==1
        % Definition of set based on percental support of density:---------
                           pik_full=pik_full/sum(pik_full);
                           % Find phenotype
                           pik_2D=reshape(pik_full,N+1,N+1);
                           [maxValue, linearIndexesOfMaxes] = max(pik_2D(:));
                           [rowsOfMaxes colsOfMaxes] = find(pik_2D == maxValue);
                             
                           % sets
                          [dens_sort,idx]=sort(pik_full,'descend');
                          last_idx=find(cumsum(dens_sort)>PercTH); 
                         
                          if rowsOfMaxes > colsOfMaxes
                            setHL=[]; 
                            setHL=idx(1:last_idx);
                            Set='setHL';
                          % save set
                            filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                            save(filename, 'setHL');
                          elseif abs(rowsOfMaxes-colsOfMaxes)<10
                            setLL=[]; 
                            setLL=idx(1:last_idx);
                            Set='setLL';
                          % save set
                            filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                            save(filename, 'setLL');
                          elseif rowsOfMaxes < colsOfMaxes
                            setLH=[];
                            setLH=idx(1:last_idx);
                            Set='setLH';
                          % save set
                            filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                            save(filename, 'setLH');
                          end
    end
    if k==2
        % Definition of set based on percental support of density:---------
                          pik_full=pik_full/sum(pik_full);
                          % Find phenotype
                          pik_2D=reshape(pik_full,N+1,N+1);
                          [maxValue, linearIndexesOfMaxes] = max(pik_2D(:));
                          [rowsOfMaxes colsOfMaxes] = find(pik_2D == maxValue);
                            
                          [dens_sort,idx]=sort(pik_full,'descend');
                          last_idx=find(cumsum(dens_sort)>PercTH);
                          
                          if rowsOfMaxes > colsOfMaxes
                          setHL=[]; 
                          setHL=idx(1:last_idx);
                          Set='setHL';
                          % save set
                          filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                          save(filename, 'setHL');
                          elseif abs(rowsOfMaxes-colsOfMaxes)<10
                          setLL=[]; 
                          setLL=idx(1:last_idx);
                          Set='setLL';
                           % save set
                          filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                          save(filename, 'setLL');
                          elseif rowsOfMaxes < colsOfMaxes
                          setLH=[];
                          setLH=idx(1:last_idx);
                          Set='setLH';
                          % save set
                            filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                            save(filename, 'setLH');
                          end
    end
    if k==3
        % Definition of set based on percental support of density:---------
                          pik_full=pik_full/sum(pik_full);
                          % Find phenotype
                          pik_2D=reshape(pik_full,N+1,N+1);                          
                          [maxValue, linearIndexesOfMaxes] = max(pik_2D(:));
                          [rowsOfMaxes colsOfMaxes] = find(pik_2D == maxValue);
                        
                          [dens_sort,idx]=sort(pik_full,'descend');
                          last_idx=find(cumsum(dens_sort)>PercTH);
                           
                        if rowsOfMaxes > colsOfMaxes
                            setHL=[]; 
                            setHL=idx(1:last_idx);
                            Set='setHL';
                        % save set
                            filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                            save(filename, 'setHL');
                        elseif abs(rowsOfMaxes-colsOfMaxes)<10
                            setLL=[]; 
                            setLL=idx(1:last_idx);
                            Set='setLL';
                           % save set
                            filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                            save(filename, 'setLL');
                         elseif rowsOfMaxes < colsOfMaxes
                            setLH=[];
                            setLH=idx(1:last_idx);
                            Set='setLH';
                            % save set
                            filename = fullfile(filepath, sprintf('SetCase%d%s.mat', Case,Set));
                            save(filename, 'setLH');
                        end
    end
end