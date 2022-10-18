%% TPT: Define the start and end sets for the TPT analysis
function []=TPTCases(Case,Q,pi,filepath,Pc,chi_crisp,tau)
%% Description:
% Function loads the sets of metastable states (i.e., phenotypes) from
% PCCA+ clustering and defines start and end sets for different transition
% paths. Finally it calculates the transition path flows between the
% metastable sets, from which percentaul paths transitions can be
% calculated (per-hand using probability tree diagrams).
%% Implemented by Anna S. Frank (anna-simone.frank@uib.no)
% load general information
Input;
myFolder = pwd;
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

switch Case
    %% Case 1: Bistability: HH and LH
    case {1}
        diary(fileTPTBi)
        Y=sprintf('Bistability -- Case %d: ', Case);
        disp(Y)
        % load data
        Filename1=sprintf('SetCase%dsetHL.mat',Case);
        Filename2=sprintf('SetCase%dsetLH.mat',Case);
        File1=fullfile(filepath,Filename1);
        File2=fullfile(filepath,Filename2);
        X=load(File1);
        Y=load(File2);
        %% ==================================================================
        %% Case 1 bistability -- Change manually set A and B to vary flow directions flow directions
           setA=X.setHL;  
           setB=Y.setLH;
        %% ==================================================================
        %% TPT prob. functions here
        % Run files for TPT analysis
        [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
        %Save file (optional)
        %filename='ComittorVars.mat';
        % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
        
        %TPT statistics
        [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probBi(qp,qm,setC,pi,Q',setA,setB,n,N,Case,filepath);
         X = sprintf('tAB = %d',tAB); disp(X)
         X = sprintf('kAB1 = %d',kAB1); disp(X)
         X = sprintf('kAB2 = %d',kAB1); disp(X)
         X = sprintf('ZAB = %d',ZAB); disp(X)
         
        %% TPT flow analysis for bistable case 1:
        prob2(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,tau)
        disp('============================================================')
        diary off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Case 2: Bistability HL and LH
         case {2}
        diary(fileTPTBi)
        Y=sprintf('Bistability -- Case %d: ', Case);
        disp(Y)
        % load data
        Filename1=sprintf('SetCase%dsetHL.mat',Case);
        Filename2=sprintf('SetCase%dsetLH.mat',Case);
        File1=fullfile(filepath,Filename1);
        File2=fullfile(filepath,Filename2);
        X=load(File1);
        Y=load(File2);
        %% ==================================================================
        %% Case 2 bistability -- Change manually set  A (START set) and B (END set) to vary flow directions
        setB=X.setHL;
        setA=Y.setLH;
        %% ==================================================================
        %% TPT prob functions here
        % Run files for TPT analysis
        [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
        %Save file (optional)
        %filename='ComittorVars.mat';
        % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
        
        % TPT statistics
        [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probBi(qp,qm,setC,pi,Q',setA,setB,n,N,Case,filepath);
         X = sprintf('tAB = %d',tAB); disp(X)
         X = sprintf('kAB1 = %d',kAB1); disp(X)
         X = sprintf('kAB2 = %d',kAB1); disp(X)
         X = sprintf('ZAB = %d',ZAB); disp(X)
         
        % TPT flow analysis for bistable case 2:
        prob2(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,tau)
        disp('============================================================')
        diary off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Case 3: Bistability LL and LH
         case {3}
        diary(fileTPTBi)
        Y=sprintf('Bistability -- Case %d: ', Case);
        disp(Y)
        % load data
        Filename1=sprintf('SetCase%dsetLL.mat',Case);
        Filename2=sprintf('SetCase%dsetLH.mat',Case);
        File1=fullfile(filepath,Filename1);
        File2=fullfile(filepath,Filename2);
        X=load(File1);
        Y=load(File2);
        %% ==================================================================
        %% Case 3 bistability -- Change manually set  A (START set) and B (END set) to vary flow directions
        setA=X.setLL;
        setB=Y.setLH;
        %% ==================================================================
        %% TPT prob functions here
        % Run files for TPT analysis
        [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
        %Save file (optional)
        %filename='ComittorVars.mat';
        % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
        
        %TPT statistics
        [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probBi(qp,qm,setC,pi,Q',setA,setB,n,N,Case,filepath);
         X = sprintf('tAB = %d',tAB); disp(X)
         X = sprintf('kAB1 = %d',kAB1); disp(X)
         X = sprintf('kAB2 = %d',kAB1); disp(X)
         X = sprintf('ZAB = %d',ZAB); disp(X)
         
        %% TPT flow analysis for bistable case 3:
        prob2(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,tau)
        disp('============================================================')
        diary off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Case 4: Bistability LL and HL
         case {4}
        diary(fileTPTBi)
        Y=sprintf('Bistability -- Case %d: ', Case);
        disp(Y)
        % load data
        Filename1=sprintf('SetCase%dsetLL.mat',Case);
        Filename2=sprintf('SetCase%dsetHL.mat',Case);
        File1=fullfile(filepath,Filename1);
        File2=fullfile(filepath,Filename2);
        X=load(File1);
        Y=load(File2);
        %% ==================================================================
        %% Case 4 bistability -- Change manually set A (START set) and B (END set) to vary flow directions
        setA=X.setLL;
        setB=Y.setHL;
        %% ==================================================================
        %% TPT prob functions here
        % Run files for TPT analysis
        [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
        %Save file (optional)
        %filename='ComittorVars.mat';
        % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
        
        % TPT Statistics
        [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probBi(qp,qm,setC,pi,Q',setA,setB,n,N,Case,filepath);
         X = sprintf('tAB = %d',tAB); disp(X)
         X = sprintf('kAB1 = %d',kAB1); disp(X)
         X = sprintf('kAB2 = %d',kAB1); disp(X)
         X = sprintf('ZAB = %d',ZAB); disp(X)
         
        % TPT flow analysis for bistable case 4:
        prob2(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,tau)
        disp('============================================================')
        diary off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Case 5: Tristable case: There are 6 subcases 
    case 5 
         % Load data
         Filename1=sprintf('SetCase%dsetLL.mat',Case);
         Filename2=sprintf('SetCase%dsetLH.mat',Case);
         Filename3=sprintf('SetCase%dsetHL.mat',Case);
         File1=fullfile(filepath,Filename1);
         File2=fullfile(filepath,Filename2);
         File3=fullfile(filepath,Filename3);
         X=load(File1);
         Y=load(File2);
         Z=load(File3);
       
         %Select which subcase to run for TPT analysis in Tristable situation      
         SubC=input('Select a subcase for tristable situation (Options 1 to 6): ')
        
         switch SubC
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Case 1: from LL to combined set LH-HL
            case 1
                diary(fileTPTTriSc1)
                disp('Tristability: SubCase 1')
                %% ==================================================================
                %% Define start and end sets for TPT: Set A: START; Set B: END
                setA=X.setLL;
                setB1=Y.setLH;
                setB2=Z.setHL;
                setB=union(setB1,setB2);
                %% ==================================================================
                % Run files for TPT analysis
                [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
                %Save file (optional)
                %filename='ComittorVars.mat';
                % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
                
                % TPT statistics
                [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probTri(qp,qm,setC,pi,Q',setA,setB,n,N,Case,SubC,filepath);
                X = sprintf('tAB = %d',tAB); disp(X)
                X = sprintf('kAB1 = %d',kAB1); disp(X)
                X = sprintf('kAB2 = %d',kAB1); disp(X)
                X = sprintf('ZAB = %d',ZAB); disp(X)
              
                % TPT flow analysis for specific subcase
                prob3(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,setB1,setB2,tau)
                diary off
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% CASE 2: from combined set LL-LH to HL (Not relevant -- I think)
            case 2
                diary(fileTPTTriSc2)
                disp('Tristability: SubCase 2')
                %% ==================================================================
                %% Define start and end sets for TPT: Set A: START; Set B: END
                setA=X.setLL;
                setB1=Z.setHL;
                setB2=Y.setLH;
                setA=union(setA,setB2);
                setB=setB1;
                %% ==================================================================
               % Run files for TPT analysis
               [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
               %Save file (optional)
               %filename='ComittorVars.mat';
               % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
               
               %% TPT statistics
               [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probTri(qp,qm,setC,pi,Q',setA,setB,n,N,Case,SubC,filepath);
               X = sprintf('tAB = %d',tAB); disp(X)
               X = sprintf('kAB1 = %d',kAB1); disp(X)
               X = sprintf('kAB2 = %d',kAB1); disp(X)
               X = sprintf('ZAB = %d',ZAB); disp(X)
              
               % Run TPT flow analysis for specific subcase
               prob3(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,setB1,setB2,tau)
               diary off
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %% CASE 3: from LL to HL, over set LH
            case 3
                diary(fileTPTTriSc3)
                disp('Tristability: SubCase 3')
                %% ==================================================================
                %% Define start and end sets for TPT: Set A: START; Set B: END
                setA=X.setLL;
                setB1=Z.setHL;
                setB2=Y.setLH;
                setB=setB1;
                %% ==================================================================
                % Run files for TPT analysis
                [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
                %Save file (optional)
                %filename='ComittorVars.mat';
                % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
                
                % TPT statistics
                [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probTri(qp,qm,setC,pi,Q',setA,setB,n,N,Case,SubC,filepath);
                X = sprintf('tAB = %d',tAB); disp(X)
                X = sprintf('kAB1 = %d',kAB1); disp(X)
                X = sprintf('kAB2 = %d',kAB1); disp(X)
                X = sprintf('ZAB = %d',ZAB); disp(X)
              
                % TPT flow analysis for specific subcase
                prob4(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,setB1,setB2,tau) %correct
                diary off
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% CASE 4: from LL to LH over set HL
            case 4
                diary(fileTPTTriSc4)
                disp('Tristability: SubCase 4')
                %% ==================================================================
                %% Define start and end sets for TPT: Set A: START; Set B: END
                setA=X.setLL;
                setB1=Y.setLH;
                setB2=Z.setHL;
                setB=setB1;
                %% ==================================================================
                % Run files for TPT analysis
                [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
                %Save file (optional)
                %filename='ComittorVars.mat';
                % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
                
                %TPT Staitsics
                [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probTri(qp,qm,setC,pi,Q',setA,setB,n,N,Case,SubC,filepath);
                X = sprintf('tAB = %d',tAB); disp(X)
                X = sprintf('kAB1 = %d',kAB1); disp(X)
                X = sprintf('kAB2 = %d',kAB1); disp(X)
                X = sprintf('ZAB = %d',ZAB); disp(X)
              
                % TPT flow analysis for specific subcase
                prob4(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,setB1,setB2,tau) % correct
                diary off
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% CASE 5: from HL to LH over set LL
            case 5
                diary(fileTPTTriSc5)
                disp('Tristability: SubCase 5')
                %% ==================================================================
                %% Define start and end sets for TPT: Set A: START; SetB: END
                setA=Z.setHL;
                setB1=Y.setLH;
                setB2=X.setLL;
                setB=setB1; % original B1
                %% ==================================================================
                % Run files for TPT analysis
                [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
                %Save file (optional)
                %filename='ComittorVars.mat';
                % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
 
                % TPT statistics
                [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probTri(qp,qm,setC,pi,Q',setA,setB,n,N,Case,SubC,filepath);
                X = sprintf('tAB = %d',tAB); disp(X)
                X = sprintf('kAB1 = %d',kAB1); disp(X)
                X = sprintf('kAB2 = %d',kAB1); disp(X)
                X = sprintf('ZAB = %d',ZAB); disp(X)
               
                % TPT flow analysis for specific subcase
                prob4(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,setB1,setB2,tau)
                diary off
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% CASE 6: from LH to HL over set LL
            case 6
                diary(fileTPTTriSc6)
                disp('Tristability: SubCase 6')
                %% ==================================================================
                %% Define start and end sets for TPT: Set A: START; Set B: END
                setA=Y.setLH;
                setB1=Z.setHL;
                setB2=X.setLL;
                setB=setB1;
                %% ==================================================================
                % Run files for TPT analysis
                [qp,qm,setC,pi,Q,setA,setB,n,N]=committor(Q',pi,setA,setB,Case,filepath);
                %Save file (optional)
                %filename='ComittorVars.mat';
                % save(filename,'Q','pi','qp','qm','setA','setB','setC','N','n','Pc','chi_crisp')
 
                % TPT statistics
                [mAB,fAB,Feff,kAB1,kAB2,ZAB,tAB]=probTri(qp,qm,setC,pi,Q',setA,setB,n,N,Case,SubC,filepath);
                X = sprintf('tAB = %d',tAB); disp(X)
                X = sprintf('kAB1 = %d',kAB1); disp(X)
                X = sprintf('kAB2 = %d',kAB1); disp(X)
                X = sprintf('ZAB = %d',ZAB); disp(X)
              
                % TPT flow analysis for specific subcase
                prob4(pi,Q',setA,setB,setC,Pc,qp,chi_crisp,fAB,setB1,setB2,tau)
                diary off
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            otherwise
                disp('No such case')
         end
end
        
