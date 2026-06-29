function [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,TrainSet] = UL_Select(l_Pop,output,Archive,RP,SQ,RQ,W,Global)


    %% Setting
    Nu = Global.N(1);
    Nl = Global.N(2);
    
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
    Au = zeros(length(A),Global.D(1));
    for i = 1:length(A)
        Au(i,:) = A{i}(1).upper_dec;
        
        if i>length(Archive)
            output = [output,A{i}(A{i}.adds==1)];
        end
    end
    
    %% Update 

     [FrontNo,~] = NDSort(output.upper_objs,[output.upper_cons,output.lower_cons],1);
      output = output(FrontNo==1);
      if length(output)>Global.N(3)
         Zmin = min(output.upper_objs);
         [Z,~]= UniformPoint(Global.N(3),2);
         [~,I] = max(1-pdist2(output.upper_objs-Zmin,Z,'cosine'));
         I = unique(I);
         output = output(I);
      end
      o_udecs = output.upper_decs;
    
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

