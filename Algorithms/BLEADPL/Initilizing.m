function [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,W,TrainSet] = Initilizing(K,Global)
% Initialize upper level (UL) population and search for their optimal lower level (LL) solutions
%
    %% Setting
    Nu = Global.N(1);
    Nl = Global.N(2);
    u_M = Global.M(1);
    
    [W,K] = UniformPoint(K,u_M);
    N   = ceil(Nu/K)*K;
    S   = N/K;   
    
    %% Initilizing    
    l_Pop = cell(1,Nu);
    RP = cell(1,Nu); 
    R(1,Nu*Nl) = INDIVIDUAL;
    TrainSet = R;
    Idx_r = zeros(length(R),1);
    Nr = 0;
    
    Pop = Global.Initialization;
    Pu = Pop.upper_decs;     
    % Lower Level (LL) multiobjective optimization run for all xu,i
    for i = 1:Nu
        
        [l_Pop{i},RP{i}] = LLsearch_DPL(Pu(i,:),Global);           %Update SPi to LL Pareto-optimum

        R(Nr+1:Nr+sum(l_Pop{i}.adds==1)) = l_Pop{i}(l_Pop{i}.adds==1);
        Idx_r(Nr+1:Nr+sum(l_Pop{i}.adds==1)) = ones(sum(l_Pop{i}.adds==1),1)*i;
        Nr = Nr +sum(l_Pop{i}.adds==1);
        
        TrainSet((i-1)*Nl+1:i*Nl) = l_Pop{i};
    end
    R(Nr+1:end) = [];
    Idx_r(Nr+1:end) = [];
    
              
      %% Update 
      [FrontNo,~] = NDSort(R.upper_objs,[R.upper_cons,R.lower_cons],1);
      output = R(FrontNo==1);
      o_udecs = output.upper_decs;
    
      [~,Loc] = ismember(o_udecs,Pu,'rows');
      I = unique(Loc);
      Archive = cell(1,length(I));
      RPs = cell(1,length(I));
      for i = 1:length(I)
          Archive{i} = l_Pop{I(i)};
          RPs{i} = RP{I(i)};
      end
      [~,Idx_oA] = ismember(Loc,I);
      
      
     [u_Pop,~,Idx_ul,l_Pop,Arch_partion] = Update_v10(R,Idx_r,l_Pop,W,S);      
     
end