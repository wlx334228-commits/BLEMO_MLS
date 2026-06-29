function [Pop,Upper_Data_Pop] = global_initialize(n_s_0,N_l)
    Global = GLOBAL.GetObj();
    U_P = Global.Initialization('upper',n_s_0);
    Pop = [];
    for i=1:n_s_0
        Data_i.upper_dec = U_P(i).upper_dec;
        Data_i.Population = Global.Evaluate(Global.Initialization(U_P(i).upper_dec,N_l));
        [~,~,Data_i.FrontNo,Data_i.CrowdDis] = NSGAII_Update(Data_i.Population,N_l,'lower');
        Pop = cat(2,Pop,Data_i);
    end
    Upper_Data_Pop = Upper_Data(Pop);
end