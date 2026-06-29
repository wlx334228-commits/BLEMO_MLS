function cG_BLEMO(Global)
% <algorithm> <N>
% multiobjective bilevel evolutionary algorithm based on multiple quadratic fibers (mf-BLEA_CGAN)

%------------------------------- Reference --------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    warning off
    Nu = Global.N(1);
    Nl = Global.N(2);
    K = 5;
    v = 1;
    
    u_PF = Global.problem.PF();
    if ~isempty(u_PF)
        tRP = 1.1*max(u_PF);
    end
    

    cgan =[];
    Record = cell(1,10);
    Record_max = zeros(10,Global.M(1));
    HVs = zeros(1,10);
    
    notermination = true;    
    gen = 1;
    
    N0 = 5;
    N = N0;

%     tmax = 100;
    tm = 1;
    %% Initialize upper level (UL) population and search for their optimal lower level (LL) solutions


    [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,W,TrainSet] = Initilize_v20(K,Global);
   
    L = cat(2,l_Pop{:});
    L = L(L.lower_adds==1);
    [FNo,~] =NDSort(L.upper_objs,[L.upper_cons,L.lower_cons],1);
    L = L(FNo==1);
    Record{1} = L;%output
    Record_max(1,:) = max(L.upper_objs);%output
    
    NFEs = Global.FEs;
    
    if ~isempty(u_PF)
        Igd  = IGD(output.upper_objs,u_PF);
        Hv   = HV(output.upper_objs,tRP);
    else
        Igd = [];
        Hv  = [];
    end
    
    %% Optimization
    while Global.NotTermination(output,true)
        %% CGAN
         if length(TrainSet)>0.7*Global.N(3)
         if isempty(cgan)
%             Notcompelted = true;
%             while Notcompelted
            try
                cgan = CGANv9(TrainSet.upper_decs,TrainSet.lower_decs,Global);
                tmax = 100;
%                 Notcompelted =false;
            catch
                tmax = 1;
            end
%             end
         elseif v<0 && tmax>1
            cgan = cgan.train(TrainSet.upper_decs,TrainSet.lower_decs,2*Global.D(2)); %5*Global.l_D
         end
         else
            tmax = 1; 
         end
        %% Generate the upper level Offspring
        
        [SQ,RQs,TrainSet] = Reproduce_v20(N,tmax,u_Pop,Idx_ul,l_Pop,Arch_partion,cgan,W,Global,TrainSet,output);
        
%% Combine parents and offspring, Update the Population
        if ~isempty(SQ)                          
            [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs] = UL_Select_m2m_v20(l_Pop,output,Archive,RPs,SQ,RQs,W,Global);

%             if mod(gen,5)
%                 W = W_update(R,Nu,Nl,K);
%             end 
        end    
        
        Global.upper_Output(u_Pop,2);
            
        %% Refinesearch
        
        [output,Idx_oA,Archive,RPs] = RefineSearch_v20(output,Idx_oA,Archive,RPs,tm,Global);
%         if sum([Global.upper_FEs,Global.lower_FEs]) >= sum(Global.maxFEs)
        if all([Global.upper_FEs,Global.lower_FEs] >= Global.maxFEs)
            Global.NotTermination(output,false);
        end

        gen = gen + 1;
        
        
        %% Termination_check
        % Record the last 10 generations of populations
        k = mod(gen,10);
        if k == 0
            k = 10;
        end
        L = cat(2,l_Pop{:});
        L = L(L.lower_adds==1);
        [FNo,~] =NDSort(L.upper_objs,[L.upper_cons,L.lower_cons],1);
        L = L(FNo==1);
        Record{k} = L;%output
        Record_max(k,:) = max(L.upper_objs);%output
         
        
        ReferencePoint = max(Record_max); 
        if sum(HVs)==0|| sum(ReferencePoint~=ReferencePoint0)>0
            ReferencePoint0 = ReferencePoint;
            for i = 1 : min(gen,10)
                HVs(i) = HV(Record{i}.upper_objs,ReferencePoint); 
            end
         else
            HVs(k) = HV(Record{k}.upper_objs,ReferencePoint);
         end
            
         HV_max = max(HVs(1 : min(gen,10)));
         HV_min = min(HVs(1 : min(gen,10)));
         
         if gen > 10    
            notermination = (HV_max - HV_min)/(HV_max + HV_min) > 0.001; 
            N =max(N,N0 + max(0,ceil(-log10(max((HV_max - HV_min)/(HV_max + HV_min),0.001)))*5));
        end
            
%          tm = max(1,ceil(-log10(max((HV_max - HV_min)/(HV_max + HV_min),0.001)))*15);
         
            index =mod([k+9,k+8],10);
            site  = index ==0;
            index(site) =10;
            
            Ir1 = abs((HVs(k)-HVs(index(1)))/(HV_max + HV_min));

            Ir2 = abs((HVs(index(1))-HVs(index(2)))/(HV_max + HV_min));
           
%             tmax = min(100,round((tmax+Ir1)^(Ir1/Ir2)));
%             if tmax ==1 && rand<=0.5
%                 tmax = 2;
%             end
            
            
%             if Ir1 < 0.001 %||Ir1 < 0.001
%                 tmax = 1;
%             else
                Ir = Ir1/Ir2;
%                 Itemp = sign(Ir-1)*min(sign(Ir-1)*Ir,sign(Ir-1)*ceil(Ir));
%                 v =v*sign(log2(Ir));
                 if Ir<1
                    v=-1*v;
                 end
%                 Itemp = Ir^(v);
                Itemp = max(0,1+v*abs(Ir-1));
                
%                 tmax = min(100,max(floor(tmax + 0.2*tmax*min(5,log2(Itemp))),1));
%                 tmax = min(100,max(floor(tmax + tmax*(5*log2(Itemp)/log2(2))),1));
                tmax = min(100,max(floor(tmax + tmax*(5^(sign(v+1))*log2(Itemp)/log2(3))),1));
%             end
        
        
%         NFEs = [NFEs;Global.u_evaluated,Global.l_evaluated];
%         if ~isempty(u_PF)
%             Igd  = [Igd;IGD(output.upper_objs,u_PF)];
%             Hv   = [Hv;HV(output.upper_objs,tRP)];
%         end
%         newP = cell(1);
%         newP{1} = output;
%         outputs = [outputs,newP];
%         
%         
%         Output{1} = NFEs;
%         Output{2} = output;
%         Output{3} = Igd;
%         Output{4} = Hv;
%         Output{5} = outputs;
        %% Display
        
%         if gen>2
%             delete(subplot(2,2,1));
%             delete(subplot(2,2,3));
%         end
%         subplot(2,2,1);       
%         title('Upper-level Obj')
%         for i = 1:Nu
%             Draw(l_Pop{i}.upper_objs);
%         end
%         Draw(u_PF,'rs');
%         I = sum(max(0,u_Pop.upper_cons),2)<=0;
%         Draw(u_Pop(I).upper_objs,'go');
%         Draw(u_Pop(~I).upper_objs,'ro');
%         
%         
%         subplot(2,2,3);        
%         title('Upper-level Dec')
%         for i = 1:Nu
%             Draw(l_Pop{i}(1).upper_dec);
%         end
%         if Global.u_D>3
%             for i = 1:Nu
%                 Draw(l_Pop{i}(1).upper_dec);
%             end
%             Draw(u_Pop(I).upper_decs,'g');
%             Draw(u_Pop(~I).upper_decs,'r');
%         else
%             for i = 1:Nu
%                 plot(l_Pop{i}(1).upper_dec,'o');
%             end
%             plot(u_Pop(I).upper_decs,'go');
%             plot(u_Pop(~I).upper_decs,'ro');
%         end
%         hold off;
%         pause(0.00001);
        
    end
end