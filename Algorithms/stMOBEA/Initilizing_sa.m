function [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,W] = Initilizing_sa(K,Global)
% Initialize upper level (UL) population and search for their optimal lower level (LL) solutions
%
    %% Setting
    Nu = Global.Nu;
    Nl = Global.Nl;
    u_M = Global.u_M;
    
    [W,K] = UniformPoint(K,u_M);
    N   = ceil(Nl/K)*K;
    S   = N/K;   
    
    %% Initilizing    
    l_Pop = cell(1,Nu);
    RP = cell(1,Nu); 
    R(1,Nu*Nl) = INDIVIDUAL;
    Idx_r = zeros(length(R),1);
    Nr = 0;

%     Pu = round(rand(Nu,Global.u_D));  
    Pu = rand(Nu,Global.u_D).*(Global.u_upper-Global.u_lower)+Global.u_lower;  
    
%     Dis = pdist2(Pu,Pu);
%     order = zeros(Nu,1);
%     order(1) = 1;  
%     neighbor = zeros(Nu,1);
%     neighbor(1) = 1;
%     for i = 2:Nu
%         index1 = find(order<1);
%         index2 = find(order>0);
%         cDis  = Dis(index1,index2);
%         [~,I] = min(cDis(:));
%         I2 = ceil(I/length(index1));
%         order(index1(I-(I2-1)*length(index1))) = i;
%         neighbor(index1(I-(I2-1)*length(index1)))=index2(I2);
%     end
    
    % Lower Level (LL) multiobjective optimization run for all xu,i
    for i = 1:Nu
%         I = find(order==i);
%         if order(I)~=neighbor(I)
%             Pt = l_Pop{neighbor(I)};
%         else
            Pt = [];
%         end
        [l_Pop{I},RP{I}] = LLsearch_st(Pu(I,:),Pt,Global);           %Update SPi to LL Pareto-optimum

        
%         l_PF = Global.problem.lower_PF(Pu(i,:));                  
%         Draw(l_PF,'rs');
%         Draw(Temp_SP{i}.lobjs);

        R(Nr+1:Nr+sum(l_Pop{I}.adds==1)) = l_Pop{I}(l_Pop{I}.adds==1);
        Idx_r(Nr+1:Nr+sum(l_Pop{I}.adds==1)) = ones(sum(l_Pop{I}.adds==1),1)*I;
        Nr = Nr +sum(l_Pop{I}.adds==1);
        
    end
    R(Nr+1:end) = [];
    Idx_r(Nr+1:end) = [];
    
              
      %% Update 
      [FrontNo,~] = NDSort(R.uobjs,[R.ucons,R.lcons],1);
      output = R(FrontNo==1);
      o_udecs = output.udecs;
    
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