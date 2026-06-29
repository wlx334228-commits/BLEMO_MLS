function [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,TrainSet] = UL_Select_st(l_Pop,output,Archive,RP,SQ,RQ,W,Global)
%UNTITLED4 此处显示有关此函数的摘要
%   Base on NSGA III

    %% Setting
    Nu = Global.Nu;
    Nl = Global.Nl;
    
    K = size(W,1);
    N   = ceil(Nu/K)*K;
    S   = N/K;
    
    %% Loading
    SR = [l_Pop,SQ];
    R(1,length(SR)*Nl) = INDIVIDUAL;
    TrainSet = R;
    Idx_r = zeros(length(R),1);
    Nr = 0;
    for i = 1:length(SR)        
        R(Nr+1:Nr+sum(SR{i}.adds==1)) = SR{i}(SR{i}.adds==1);
        Idx_r(Nr+1:Nr+sum(SR{i}.adds==1)) = ones(sum(SR{i}.adds==1),1)*i;
        Nr = Nr +sum(SR{i}.adds==1);
                        
        TrainSet((i-1)*Nl+1:i*Nl) = SR{i};
    end
    R(Nr+1:end) = [];
    Idx_r(Nr+1:end) = [];
    
    A = [Archive,SQ];
    U =[RP,RQ];
    Au = zeros(length(A),Global.u_D);
    for i = 1:length(A)
        Au(i,:) = A{i}(1).udec;
        
        if i>length(Archive)
            output = [output,A{i}(A{i}.adds==1)];
        end
    end
    
    %% Update 

     [FrontNo,~] = NDSort(output.uobjs,[output.ucons,output.lcons],1);
      output = output(FrontNo==1);
      if length(output)>Nu*Nl
         Zmin = min(output.uobjs);
         [Z,~]= UniformPoint(Nu*Nl,2);
         [~,I] = max(1-pdist2(output.uobjs-Zmin,Z,'cosine'));
         output = output(I);
      end
      o_udecs = output.udecs;
    
      [~,Loc] = ismember(o_udecs,Au,'rows');
      I = unique(Loc);
      Archive = cell(1,length(I));
      RPs = cell(1,length(I));
      for i = 1:length(I)
          Archive{i} = A{I(i)};
          RPs{i} = U{I(i)};
      end
      [~,Idx_oA] = ismember(Loc,I);
    
    
    [u_Pop,~,Idx_ul,l_Pop,Arch_partion] = Update_v10(R,Idx_r,SR,W,S);  
   
end

