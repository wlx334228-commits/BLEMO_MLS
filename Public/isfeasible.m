function [feasible,CV] = isfeasible(Pop,type,epsilion)
    feasible(length(Pop),1) = false;
     if nargin > 2
         switch type
             case 'upper'
                 CV = sum(max(0,Pop.upper_cons - epsilion),2);
             case 'lower'
                 CV = sum(max(0,Pop.lower_cons - epsilion),2);
         end
     elseif nargin > 1
         switch type
             case 'upper'
                 CV = sum(max(0,Pop.upper_cons),2);
             case 'lower'
                 CV = sum(max(0,Pop.lower_cons),2);
         end
%      else
%          CV = sum(max(0,[Pop.upper_cons,Pop.lower_cons]),2);
     end
     feasible(CV == 0) = true;
end