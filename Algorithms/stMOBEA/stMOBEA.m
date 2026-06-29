function stMOBEA(Global)
% <algorithm> <N>
% SABLEA-PM

%------------------------------- Reference --------------------------------

%------------------------------- Copyright --------------------------------

%--------------------------------------------------------------------------

    %% Parameter setting
    warning off
    
    Nu = Global.N(1);
    Nl = Global.N(2);
    K = 5;
    
%     if Global.problem.Parameter.area<6000
%         uFEs_max = 2500;
%     elseif Global.problem.Parameter.area<8000
%         uFEs_max = 500;
%     else
%         uFEs_max = 700;
%     end
    
    Record = cell(1,10);
    Record_max = zeros(10,Global.M(1));
    hv = zeros(1,10);
    
    notermination = true;    
    gen = 1;
    
    u_PF = Global.problem.PF();
    if ~isempty(u_PF)
        tRP = 1.1*max(u_PF);
    end
    
    N0 = 5;
    N = N0;

    tmax = 40;
    tm = 1;
    %% Initialize upper level (UL) population and search for their optimal lower level (LL) solutions


    [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,W] = Initilizing_st(K,Global);
    TrainSet = l_Pop;
   
    %% Optimization
    while Global.NotTermination(output,true)
        %% train model
        planemodel = Plane_model(TrainSet);
        

        %% Generate the upper level Offspring
        
        [SQ,RQs] = Reproduce_st(tmax,u_Pop,Idx_ul,l_Pop,Arch_partion,TrainSet,planemodel,W,Global);
        TrainSet = [l_Pop,SQ];
        
        %% Combine parents and offspring, Update the Population
        if ~isempty(SQ)                          
            [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,~] = UL_Select_BSP(l_Pop,output,Archive,RPs,SQ,RQs,W,Global);
        end    
            
        %% Refinesearch
        
        [output,Idx_oA,Archive,RPs] = RefineSearch_st(output,Idx_oA,Archive,RPs,tm,Global);

        if all([Global.upper_FEs,Global.lower_FEs] >= Global.maxFEs)
            Global.NotTermination(output,false);
        end
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
% %         Draw(u_PF,'rs');
%         I = sum(u_Pop.upper_cons,2)<=0;
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