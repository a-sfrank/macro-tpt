function [r,s]=coords(node,N) 
%computes 2D coordinates of a node with index node
%% Implemented by Susanna Röblitz (susanna.roblitz@uib.no)
    r=floor((node-1)/N);
    s=node-1-r*(N);

end
   