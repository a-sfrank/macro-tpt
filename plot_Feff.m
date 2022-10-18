
function plot_Feff(Case,SubC)
% Description: Function plots the transition flows between phenotype states
% in the tristable situation
%% Run this function separately after the main_pcca_macrophage.m file has run.
%% Implemented by Anna S. Frank (anna-simone.frank@uib.no) and Susanna Röblitz (susanna.roblitz@uib.no)
filepath='.\Results_Main_PCCA'; % put in the correct filepath, check that it works

    %Load data
    fname=sprintf('case%d_SubCase%d_Feff.mat',Case,SubC);
    x=load(fname)
    Feff=x.Feff;
    setA=x.setA;
    setB=x.setB;
    
    %Plot data
    [m,n]=size(Feff);
    N=round(sqrt(m))
    
    sf=max(max(Feff)); %Scale factor
    %color gradient: a*Feff+b
    a=0.5/(eps-sf);
    b=sf*0.5/(sf-eps);
    
    
    [ii,jj]=find(abs(Feff)>eps);
    length(setA)
   length( 1:length(setA))
   
    figure(33)
    for k=1:length(setA)
        [i,j]=coords(setA(k),101);
        plot(i,j,'r*');
        hold on
    end
    for k=1:length(setB)
        [i,j]=coords(setB(k),101);
        plot(i,j,'bx');
        hold on
    end
    
    for k=1:length(ii)
        k
        % Plotbefehl fur Kanten
        r=floor((ii(k)-1)/(N));
        s=ii(k)-1-r*(N);
        t=floor((jj(k)-1)/(N));
        u=jj(k)-1-t*(N);
        z=max(0,a*Feff(ii(k),jj(k))+b);
        h=quiver(r,s,t-r,u-s,0,'k','Color',[z z z]); 
        xlim([0 100]); ylim([0 100])
        hold on
    end
    axis equal
    xlim([0 100]); ylim([0 100])
    filepath='.\Results_Main_PCCA';
    fname = sprintf('Fig_33_Case_%d_SubC_%d.fig',Case,SubC);
    saveas(h, fullfile(filepath,fname));
    xlabel('x_1')
    ylabel('x_2')
    
    fname2 = sprintf('Fig_33_Case_%d_SubC_%d.eps',Case,SubC);
    saveas(h, fullfile(filepath,fname2));
    xlabel('x_1')
    ylabel('x_2')
    
    fname3 = sprintf('Fig_33_Case_%d_SubC_%d.png',Case,SubC);
    saveas(h, fullfile(filepath,fname3));
    xlabel('x_1')
    ylabel('x_2')
end