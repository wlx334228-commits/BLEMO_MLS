function [I_epsilon,rate] = Cal_I_epsilon(PopObj1,PopObj2)

    N1 = size(PopObj1,1);
    N2 = size(PopObj2,1);
    Min_obj = min([PopObj1;PopObj2],[],1);
    
    PopObj1 = PopObj1 - Min_obj;
    PopObj2 = PopObj2 - Min_obj;
    
    I_epsilon = zeros(1,N1);
    rate = zeros(1,N1);
    for i=1:N1
%         I_epsilon(i) = min(max(PopObj2./max(repmat(PopObj1(i,:),N2,1),eps),[],2));
        I_epsilon(i) = min(max(repmat(PopObj1(i,:),N2,1)./max(PopObj2,eps),[],2));
        rate(i) = sum(max(repmat(PopObj1(i,:),N2,1)./max(PopObj2,eps),[],2)<1);
    end
    
end