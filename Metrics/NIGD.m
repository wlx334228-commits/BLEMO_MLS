function Score = NIGD(PopObj,PF)
% <metric> <min>
% Inverted generational distance

%------------------------------- Reference --------------------------------
% C. A. Coello Coello and N. C. Cortes, Solving multiobjective optimization
% problems using an artificial immune system, Genetic Programming and
% Evolvable Machines, 2005, 6(2): 163-190.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    M = size(PopObj,1);N=size(PF,1);
    PopObj = (PopObj-repmat(min(PF),M,1))./(repmat(max(PF)-min(PF),M,1));
    PF=(PF-repmat(min(PF),N,1))./(repmat(max(PF)-min(PF),N,1));
%     Interc=max(PF,[],1);
%     PF=PF./repmat(Interc,size(PF,1),1);
%     PopObj=PopObj./repmat(Interc,size(PopObj,1),1);
    Distance = min(pdist2(PF,PopObj),[],2);
    Score    = mean(Distance);
end