function Score = IGD_epsilon(PopObj,PF)
% <metric> <min>
% Inverted generational distance

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group
     M = size(PopObj,1); N=size(PF,1);I = zeros(N,1);
      PopObj = (PopObj-repmat(min(PF),M,1))./(repmat(max(PF)-min(PF),M,1)+1e-4);
      PF=(PF-repmat(min(PF),N,1))./(repmat(max(PF)-min(PF),N,1)+1e-4);
     for i=1:N
            Xi = PF(i,:);
            CV = PopObj - repmat(Xi,M,1);
            CV1 = CV;
            CV1(CV1<0) = 0;
            %Fitness = max(CV1,[],2);
            Fitness = sum(CV1,2);
%            MaxCV=max(CV1,[],2);
%            DomInds=find(MaxCV<=0);
%            CV2 = CV(DomInds,:);
%           Fitness(DomInds) = -(sum(abs(CV2),2));
           I(i)=min(Fitness);
     end 
     Score = mean(I);
end