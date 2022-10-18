function [Q,Case]=assembleQ(N)
% Code assembles the CME matrix and loads different cases for parameter
% values
%% Implemented by Susanna Röblitz (susanna.roblitz@uib.no), Anna S. Frank (anna-simone.frank@uib.no)
 %% Load parameters 
    Case = input('Which parameter case do you want to run. Choose between numbers 1 to 5: ');
    [a1,a2,b1,b2,k1,k2,p1,p2,q1,q2,S1,S2,n1,n2,l1,l2]=parameterCase(Case) ;  

%% Assembly CME matrix
    Q=sparse((N+1)^2,(N+1)^2);
    for i=1:(N-1)
        for j=1:(N-1)
            stateIndex=i*(N+1)+j+1;
            state_im1=(i-1)*(N+1)+j+1;
            state_ip1=(i+1)*(N+1)+j+1;
            state_jp1=i*(N+1)+j+2;
            state_jm1=i*(N+1)+j;
            alpha1p=(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
            alpha2p=a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
            alpha1m=(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
            alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
            Q(stateIndex,state_im1)=alpha1p;
            Q(stateIndex,state_ip1)=q1*(i+1);
            Q(stateIndex,state_jm1)=alpha2p;
            Q(stateIndex,state_jp1)=q2*(j+1);
            Q(stateIndex,stateIndex)=-alpha1m-q1*i-alpha2m-q2*j;
        end
    end
    %deal with 1st and last row and column seperately
    i=0;
    for j=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_ip1=(i+1)*(N+1)+j+1;
        state_jm1=i*(N+1)+j;
        state_jp1=i*(N+1)+j+2;
        alpha2p=a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        alpha1m=(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_ip1)=q1*(i+1);
        Q(stateIndex,state_jm1)=alpha2p;
        Q(stateIndex,state_jp1)=q2*(j+1);
        Q(stateIndex,stateIndex)=-alpha1m-alpha2m-q2*j;
    end
    j=0;
    for i=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_im1=(i-1)*(N+1)+j+1;
        state_ip1=(i+1)*(N+1)+j+1;
        state_jp1=i*(N+1)+j+2;
        alpha1p=(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha1m=(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_im1)=alpha1p;
        Q(stateIndex,state_ip1)=q1*(i+1);
        Q(stateIndex,state_jp1)=q2*(j+1);
        Q(stateIndex,stateIndex)=-alpha1m-q1*i-alpha2m;
    end
    i=N;
    for j=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_im1=(i-1)*(N+1)+j+1;
        state_jm1=i*(N+1)+j;
        state_jp1=i*(N+1)+j+2;
        alpha1p=(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2p=a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        alpha1m=(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_im1)=alpha1p;
        Q(stateIndex,state_jm1)=alpha2p;
        Q(stateIndex,state_jp1)=q2*(j+1);
        Q(stateIndex,stateIndex)=-q1*i-alpha2m-q2*j;    %-alpha1m
    end
    j=N;
    for i=1:(N-1)
        stateIndex=i*(N+1)+j+1;
        state_im1=(i-1)*(N+1)+j+1;
        state_ip1=(i+1)*(N+1)+j+1;
        state_jm1=i*(N+1)+j;
        alpha1p=(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2p=a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        alpha1m=(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
        alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
        Q(stateIndex,state_im1)=alpha1p;
        Q(stateIndex,state_ip1)=q1*(i+1);
        Q(stateIndex,state_jm1)=alpha2p;
        Q(stateIndex,stateIndex)=-alpha1m-q1*i-q2*j;    %-alpha2m
    end
    i=0;
    j=0;
    stateIndex=i*(N+1)+j+1;
    state_ip1=(i+1)*(N+1)+j+1;
    state_jp1=i*(N+1)+j+2;
    alpha1m=(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_ip1)=q1*(i+1);
    Q(stateIndex,state_jp1)=q2*(j+1);
    Q(stateIndex,stateIndex)=-alpha1m-alpha2m;
    %%%%%%%%%%%%%%%%%%%%%%%
    i=0;
    j=N;
    stateIndex=i*(N+1)+j+1;
    state_ip1=(i+1)*(N+1)+j+1;
    state_jm1=i*(N+1)+j;
    alpha2p=a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    alpha1m=(a1*i^n1/(k1^n1+i^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_ip1)=q1*(i+1);
    Q(stateIndex,state_jm1)=alpha2p;
    Q(stateIndex,stateIndex)=-alpha1m-q2*j; %-alpha2m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=N;
    j=0;
    stateIndex=i*(N+1)+j+1;
    state_im1=(i-1)*(N+1)+j+1;
    state_jp1=i*(N+1)+j+2;
    alpha1p=(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2m=a2*(j^n2/(k2^n2+j^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_im1)=alpha1p;
    Q(stateIndex,state_jp1)=q2*(j+1);
    Q(stateIndex,stateIndex)=-alpha2m-q1*i;     %-alpha1m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=N;
    j=N;
    stateIndex=i*(N+1)+j+1;
    state_im1=(i-1)*(N+1)+j+1;
    state_jm1=i*(N+1)+j;
    alpha1p=(a1*(i-1)^n1/(k1^n1+(i-1)^n1)+S1)*(1./(1+(j/p2)^l2))+b1;
    alpha2p=a2*((j-1)^n2/(k2^n2+(j-1)^n2))+S2*(1/(1+(i/p1)^l1))+b2;
    Q(stateIndex,state_im1)=alpha1p;
    Q(stateIndex,state_jm1)=alpha2p;
    Q(stateIndex,stateIndex)=-q1*i-q2*j;    %-alpha1m-alpha2m
end