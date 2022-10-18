% Function contains different parameter cases for the results presented in the article
% "Macrophage phenotype transitions in a stochastic gene-regulatory netwrok model" by Anna S. Frank,
% Kamila Larripa, Hwayeon Ryu and Susanna Röblitz:
% There are four cases of bistability and one case of tristability
%% Implemented by Anna.S.Frank (anna-simone.frank@uib.no)
function [a1,a2,b1,b2,k1,k2,p1,p2,q1,q2,S1,S2,n1,n2,l1,l2]=parameterCase(Case)

switch Case
    case 1 
        disp('Case 1: Bistability -- High/High and Low/High')
       %parameter values
       P1 = 0.5;
       P2 = 1;
       A1 = 5;
       A2 = 5; 
       K1 = 1;
       K2 = 1;
       B1 = 15;
       B2 = 0.05;
       n1 = 6; 
       n2 = 6;
       l1 = 1;
       l2 = 1;
       q1 = 5;
       q2 = 1;
       s1 = 2;
       s2 = 2/12.9;     
       V = 5e-12;    %Volume of macrophage in L
       u = 0.43*1e-11;  
       %--------------------------------------------------------
       % Parameter transformation
       nA = 6e+23; %Avrogado constant
       uVnA = u*V*nA; %3e+2; conversion factor for a,b,k,p,S
       a1 = A1*uVnA;
       a2 = A2*uVnA;
       b1 = B1*uVnA;
       b2 = B2*uVnA;
       k1 = K1*uVnA;
       k2 = K2*uVnA;
       p1 = P1*uVnA;
       p2 = P2*uVnA;
      S1 = s1*uVnA;
      S2 = s2*uVnA; 
       
      case 2 
        disp('Case 2: Bistability -- Low/High and High/Low')
        %parameter values
        P1 = 0.5;
        P2 = 1;
        A1 = 15;
        A2 = 5; 
        K1 = 1;
        K2 = 1;
        B1 = 0.05;
        B2 = 0.05;
        n1 = 6; 
        n2 = 6;
        l1 = 1;
        l2 = 1;
        q1 = 5;
        q2 = 1;
        s1 = 2;
        s2 = 2/12.9;    
        V = 5e-12;            %Volume of macrophage in L
        u = 0.43*1e-11; 
        %--------------------------------------------------------
        % Parameter transformation
        nA = 6e+23; %Avrogado constant
        uVnA = u*V*nA; %3e+2; conversion factor for a,b,k,p,S
        a1 = A1*uVnA;
        a2 = A2*uVnA;
        b1 = B1*uVnA;
        b2 = B2*uVnA;
        k1 = K1*uVnA;
        k2 = K2*uVnA;
        p1 = P1*uVnA;
        p2 = P2*uVnA;
       S1 = s1*uVnA;
       S2 = s2*uVnA; 

    case 3 
        disp('Case 3: Bistability -- Low/High and Low/Low')
        % parameter values
        P1 = 0.5;
        P2 = 1;
        A1 = 15;
        A2 = 5; 
        K1 = 1;
        K2 = 1;
        B1 = 0.05;
        B2 = 0.05;
        n1 = 6; 
        n2 = 6;
        l1 = 1;
        l2 = 1;
        q1 = 5;
        q2 = 5;
        s1 = 1;
        s2 = 1/12.9;     
        V = 5e-12;    %Volume of macrophage in L
        u = 0.43*1e-11;  
        %--------------------------------------------------------
        % Parameter transformation
        nA = 6e+23; %Avrogado constant
        uVnA = u*V*nA; %3e+2; conversion factor for a,b,k,p,S
        a1 = A1*uVnA;
        a2 = A2*uVnA;
        b1 = B1*uVnA;
        b2 = B2*uVnA;
        k1 = K1*uVnA;
        k2 = K2*uVnA;
        p1 = P1*uVnA;
        p2 = P2*uVnA;
        S1 = s1*uVnA;
        S2 = s2*uVnA; 
    
     case 4 
        disp('Case 4: Bistability -- High/Low and Low/Low')
        % parameter values
        P1 = 0.5;
        P2 = 1;
        A1 = 5;
        A2 = 5; 
        K1 = 1;
        K2 = 1;
        B1 = 0.001;
        B2 = 0.05;
        n1 = 6; 
        n2 = 6;
        l1 = 1;
        l2 = 1;
        q1 = 5;
        q2 = 1;
        s1 = 2;
        s2 = 2/12.9;     
        V = 5e-12;    %Volume of macrophage in L
        u = 0.43*1e-11;  
        %--------------------------------------------------------
        % Parameter transformation
        nA = 6e+23; %Avrogado constant
        uVnA = u*V*nA; %3e+2; conversion factor for a,b,k,p,S
        a1 = A1*uVnA;
        a2 = A2*uVnA;
        b1 = B1*uVnA;
        b2 = B2*uVnA;
        k1 = K1*uVnA;
        k2 = K2*uVnA;
        p1 = P1*uVnA;
        p2 = P2*uVnA;
        S1 = s1*uVnA;
        S2 = s2*uVnA; 
        
     case 5 
        disp('Case 5: Tristability -- Low/High and High/Low and Low/Low')
        % parameter values
        P1 = 0.5;
        P2 = 1;
        A1 = 15;
        A2 = 5; 
        K1 = 1;
        K2 = 1;
        B1 = 0.05;
        B2 = 0.05;
        n1 = 6; 
        n2 = 6;
        l1 = 1;
        l2 = 1;
        q1 = 5;
        q2 = 1;
        s1 = 4;
        s2 = 4/12.9;     
        V = 5e-12;    %Volume of macrophage in L
        u = 0.43*1e-11;  
        %--------------------------------------------------------
        % Parameter transformation
        nA = 6e+23; %Avrogado constant
        uVnA = u*V*nA; %3e+2; conversion factor for a,b,k,p,S
        a1 = A1*uVnA;
        a2 = A2*uVnA;
        b1 = B1*uVnA;
        b2 = B2*uVnA;
        k1 = K1*uVnA;
        k2 = K2*uVnA;
        p1 = P1*uVnA;
        p2 = P2*uVnA;
        S1 = s1;
        S2 = s2; 
    otherwise
        disp('No such case')
end