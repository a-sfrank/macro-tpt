
function plot_Feff_bistability(Case)
% Example run: plot_Feff_bistability(1)
% Description: Function plots the transition flows between phenotype states
% in the bistable situation
%% Run this function separately after the main_pcca_macrophage.m file has run.
%% Implemented by Anna S. Frank (anna-simone.frank@uib.no) and Susanna Röblitz (susanna.roblitz@uib.no)
filepath='.\Results_Main_PCCA'; % put in the correct filepath, check that it works

    %Load data
    fname=sprintf('case%d_Feff.mat',Case);
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
   
    fname3 = sprintf('Fig_33_Case_%d.png',Case);
    title('Transition path flow directions')
    saveas(h, fullfile(filepath,fname3));
    xlabel('x_1')
    ylabel('x_2')
end