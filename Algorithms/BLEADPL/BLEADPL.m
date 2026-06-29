function BLEADPL(Global)

    K = 5;
    
    [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,W,~] = Initilizing(K,Global);
    
    tao = 10;
    Epsilon_u = 1e-2;
    
    
    Objs = {};
    
    while Global.NotTermination(output,true)
        
        [SQ,RQs] = Reproduce(u_Pop,Idx_ul,l_Pop,Arch_partion,Global);
        
        if ~isempty(SQ)                          
            [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,~] = UL_Select(l_Pop,output,Archive,RPs,SQ,RQs,W,Global);

        end 
        
        [output,Idx_oA,Archive,RPs] = RefineSearch(output,Idx_oA,Archive,RPs,Global);
        
        Global.upper_Output(u_Pop,2);
        
%         if all([Global.upper_FEs,Global.lower_FEs] > Global.maxFEs)
%             Global.NotTermination(Elite_Population,false);
%         end
        
        if ~isempty(output)
            Objs = cat(2,Objs,{output.upper_objs});
        end
        
        if length(Objs)==tao
            RefPoint = max(cat(1,Objs{:}),[],1);
            HV = cell2mat(cellfun(@(x) CalHV(x,RefPoint), Objs, 'UniformOutput',false));
            H_u = (max(HV)-min(HV))/(max(HV)+min(HV));
            if H_u <= Epsilon_u && all([Global.upper_FEs,Global.lower_FEs] > Global.maxFEs)
                Global.NotTermination(output,false);
            end
            Objs(1) = [];
        end
    end
end
