classdef MEC < PROBLEM
    methods
        function obj = MEC()
            obj.Global.M = [2 2];
            if isempty(obj.Global.D)
                N = 20;
                obj.Global.D = [3*N,N];
            else
                N = obj.Global.D(1)/3;
            end
            
            file = fullfile('Problems','MEC');
            obj.Parameter.UE.N        =   N;                                                      %the number of devices
            obj.Parameter.UE.fd       =   load([file,'\Device data\fd_',num2str(N),'.dat']);            %computational capability of the mobile device
            obj.Parameter.U.H         =   4.11*(3/(9.15*4*pi))^2*(1./(0.05*load([file,'\Device data\H_',num2str(N),'.dat']).^(-1/4))).^2;             %channel gain
%             obj.Parameter.U.H         =   load([file,'\Device data\H_',num2str(N),'.dat']);             %channel gain
            obj.Parameter.U.F         =   load([file,'\Task data\F_',num2str(N),'.dat']);               %total number of the CPU cycles
            obj.Parameter.U.D         =   load([file,'\Task data\D_',num2str(N),'.dat']);               %the size of the input data
            obj.Parameter.U.B         =   load([file,'\Task data\B_',num2str(N),'.dat']);               %the size of the input data
            obj.Parameter.UE.fc       =   2e11;                                                   %the computational capability of the cloud server
            obj.Parameter.UE.f        =   [obj.Parameter.UE.fc;obj.Parameter.UE.fd];                                         %computational capability
            obj.Parameter.UE.pt       =   1.3;                                                    %transmission power
            obj.Parameter.UE.pr       =   0.8;                                                    %reception power
            obj.Parameter.UE.Tmax     =   1.5;                                                      %maximum time constraints
            obj.Parameter.UE.Th       =   0.5;
            obj.Parameter.UE.u        =   0.8;
            parameter.W =   20*10^6;                                                %channel bandwidth
            parameter.wn=   10^(-10);                                               %background interference power
            parameter.alpha =   10^(-27);
            parameter.k  =   10^(-10);
            
            %Algorithm parameters
%             parameter.pheromoneIni      =   load([file,'\Task data\pheromoneIni_',num2str(N),'.dat']);
%             parameter.pheromoneIni      =   1/(obj.Parameter.UE.N*parameter.pheromoneIni);
%             parameter.pheromone         =   parameter.pheromoneIni*ones(obj.Parameter.UE.N,obj.Parameter.UE.N+1);
%             parameter.beta              =   2;
%             parameter.q                 =   0.9;
%             parameter.phi               =   0.1;
%             parameter.rho               =   0.1;
            
            obj.Parameter.Rmax    =   parameter.W*log2(1+obj.Parameter.UE.pt*obj.Parameter.U.H/parameter.wn);                                     %the maximum uplink data rate
            obj.Parameter.T_t     =   repmat(obj.Parameter.U.D,1,obj.Parameter.UE.N+1)./obj.Parameter.Rmax + repmat(obj.Parameter.U.B,1,obj.Parameter.UE.N+1)./obj.Parameter.Rmax;                        %the maximum transmission time
            obj.Parameter.fmin    =   repmat(obj.Parameter.U.F,1,obj.Parameter.UE.N+1)./(repmat(obj.Parameter.UE.Tmax,1,obj.Parameter.UE.N+1)-obj.Parameter.T_t-obj.Parameter.UE.Th);                           %the minimum computational resources
            obj.Parameter.M       =   (obj.Parameter.fmin>=0&obj.Parameter.fmin<=repmat(obj.Parameter.UE.f',obj.Parameter.UE.N,1));                                            %the feasible candidate execution mode
            
            obj.Parameter.Etmin = (obj.Parameter.UE.pt.*repmat(obj.Parameter.U.D,1,obj.Parameter.UE.N)+obj.Parameter.UE.pr.*repmat(obj.Parameter.U.B,1,obj.Parameter.UE.N))./obj.Parameter.Rmax(:,2:end)+(obj.Parameter.UE.pr.*repmat(obj.Parameter.U.D,1,obj.Parameter.UE.N)+obj.Parameter.UE.pt.*repmat(obj.Parameter.U.B,1,obj.Parameter.UE.N))./obj.Parameter.Rmax(:,2:end); %the minimum transmission energy consumption
            obj.Parameter.Ecmin   =   parameter.alpha.*obj.Parameter.fmin(:,2:end).^2.*repmat(obj.Parameter.U.F,1,obj.Parameter.UE.N);
            
            obj.Global.Parameter = parameter;
            obj.Global.upper_domain = [zeros(1,3*N);N*ones(1,N),20*ones(1,2*N)];
            
%             obj.Global.upper_domain = obj.Parameter.M;
%             obj.Global.lower_domain = obj.Parameter.fmin;

        end
        
        
        function s = getUpperDecs(obj,N)
            
            s = Decs(obj,unifrnd(repmat(obj.Global.upper_domain(1,:),N,1),repmat(obj.Global.upper_domain(2,:),N,1)),'upper');
            return
            M = obj.Parameter.M;
            UE = obj.Parameter.UE;
            U = obj.Parameter.U;
            Etmin = obj.Parameter.Etmin;
            Ecmin = obj.Parameter.Ecmin;
            parameter = obj.Global.Parameter;
            
            s = -1*ones(UE.N,N);
%             s = unifrnd(repmat(obj.Global.upper_domain(1,:),N,1),repmat(obj.Global.upper_domain(2,:),N,1));
            
            for k=1:N
                table = obj.Parameter.M;
            
                Rtemp0 = [];
                Rtemp1 = [];
                
                [~,selectionOrder] = sort(sum(M,2));
                for i = 1:UE.N
                    currentTask = selectionOrder(i);
                    if s(currentTask,k)==-1
                        allowMode = find(table(currentTask,:)==1);
                        if ismember(1,allowMode)
                            cloudTask = find(s(:,k) == 0);
                            archive = [cloudTask;currentTask];
                            Hsum = repmat(sum(U.H(archive,1)),size(archive,1),1)-U.H(archive,1);
                            Rtemp = parameter.W*log2(1+UE.pt*U.H(archive,1)./(parameter.wn+Hsum));
                            ftemp = U.F(archive)./(UE.Tmax-(U.D(archive)+U.B(archive))./Rtemp);
                            if sum(ftemp)>UE.f(1) || sum(ftemp<0)>0
                                allowMode(1)=[];
%                             else
%                                 Rtemp1 = Rtemp;
                            end
                            
                            existMode = s(s(:,k) > 0,k)+1;
                            for j=1:length(existMode)
                                ind = allowMode==existMode(j);
                                allowMode(ind) = [];
                            end
                            
                        end
                        
                        if ~isempty(allowMode)
                            s(currentTask,k) = randsample(allowMode,1)-1;
                        else
                            allowMode = find(table(currentTask,:)==1);
                            s(currentTask,k) = randsample(allowMode,1)-1;
                        end
                        
%                         
%                         if ~isempty(allowMode)
%                             detalE = zeros(1,size(allowMode,2));
%                             if ismember(1,allowMode)
%                                 detalE(1) = sum((UE.pt.*U.D(archive)+UE.pt.*U.B(archive))./Rtemp1)-sum((UE.pt.*U.D(cloudTask)+UE.pt.*U.B(cloudTask))./Rtemp0);
%                             else
%                                 detalE(1) = Etmin(currentTask,allowMode(1)-1)+Ecmin(currentTask,allowMode(1)-1);
%                             end
%                             detalE(2:end) = Etmin(currentTask,allowMode(2:end)-1)+Ecmin(currentTask,allowMode(2:end)-1);
%                             Eta = 1./detalE;
%                             p = parameter.pheromone(currentTask,allowMode).*Eta.^parameter.beta;
%                             if rand<parameter.q
%                                 index = find(p==max(p));
%                             else
%                                 p = p/sum(p);
%                                 pc = cumsum(p);
%                                 index = find(pc > rand);
%                             end
%                             targetMode = allowMode(index(1));
%                             s(currentTask,k) = targetMode(1)-1;
%                             if targetMode(1)~=1
%                                 table(:,targetMode(1)) = 0;
%                             else
%                                 Rtemp0 = Rtemp;
%                                 cloudTask = [cloudTask;currentTask];
%                             end
%                             % Local pheromone management
%                             obj.Global.Parameter.pheromone(currentTask,s(currentTask,k)+1) = (1-parameter.rho)*parameter.pheromone(currentTask,s(currentTask,k)+1)+parameter.rho*parameter.pheromoneIni;
%                         end
                    end
                end
            end
        end
        
        function fmin = getLowerDecs(obj,s)
            s = s';
            UE = obj.Parameter.UE;
            U = obj.Parameter.U;
            parameter = obj.Global.Parameter;
            
            N = size(s,2);
            fmin = zeros(UE.N,N);                  % computing resource
            
            for k=1:N
                %Resource allocation
                for i = 1:UE.N
                    fmin(i,k) = obj.Parameter.fmin(i,s(i,k)+1);
                end
            end
            
            for k=1:N
                cloudTask = find(s(1:UE.N,k) == 0);
                Hsum = repmat(sum(U.H(cloudTask,1)),size(cloudTask,1),1)-U.H(cloudTask,1);
                Rtemp = parameter.W*log2(1+UE.pt*U.H(cloudTask,1)./(parameter.wn+Hsum));
                fmin(cloudTask,k) = U.F(cloudTask)./(UE.Tmax-(U.D(cloudTask)+U.B(cloudTask))./Rtemp);
            end
            
            if N==1
                obj.Global.lower_domain = [];
                obj.Global.lower_domain(1,:)=fmin';
                obj.Global.lower_domain(2,:)=UE.f(s(1:UE.N)+1)';
%                 obj.Global.lower_domain(1,:)=[fmin',50*ones(size(fmin',1),2)];
%                 obj.Global.lower_domain(2,:)=[UE.f(s'+1)',100*ones(size(fmin',1),2)];
            else
                obj.Global.lower_domain = [];
%                 obj.Global.lower_domain ={[fmin',50*ones(size(fmin',1),2)],[UE.f(s'+1),100*ones(size(fmin',1),2)]};
                obj.Global.lower_domain ={fmin',UE.f(s(1:UE.N,:)'+1)};
            end
        end
        
        function Decs = Decs(obj,Decs,type)
            UE = obj.Parameter.UE;
            U = obj.Parameter.U;
            parameter = obj.Global.Parameter;
            switch type
                case 'upper'
%                     Decs = floor(Decs);
                    Decs(:,1:UE.N) = round(Decs(:,1:UE.N));
                    domain = [];
                    for i=1:UE.N
                        domain = cat(1,domain,{find(obj.Parameter.M(i,:))-1});
                    end
                    [~,selectionOrder] = sort(sum(obj.Parameter.M,2));
                    for i=1:UE.N
                        ii = selectionOrder(i);
                        currentMode = Decs(:,ii);
                        allowMode_i = domain{ii};
                        inds = ~ismember(currentMode,allowMode_i);
                        for j=1:size(Decs,1)
                            
                            allowMode_i = domain{ii};
                            
                            assigned_t = selectionOrder(1:i-1);
                            cloudTask = assigned_t(Decs(j,assigned_t) == 0);
                            archive = [cloudTask;ii];
                            
                            Hsum = repmat(sum(U.H(archive,1)),size(archive,1),1)-U.H(archive,1);
                            Rtemp = parameter.W*log2(1+UE.pt*U.H(archive,1)./(parameter.wn+Hsum));
                            ftemp = U.F(archive)./(UE.Tmax-(U.D(archive)+U.B(archive))./Rtemp);
                            if sum(ftemp)>UE.f(1) || sum(ftemp<0)>0
                                allowMode_i(1)=[];
%                             else
%                                 Rtemp1 = Rtemp;
                            end
                            
                            inds_j = ismember(allowMode_i,Decs(j,assigned_t));
                            inds_j(allowMode_i==0)=false;
                            if inds(j) || (currentMode(j)>0 && ismember(currentMode(j),Decs(j,selectionOrder(1:i-1)))) || ~ismember(Decs(j,ii),allowMode_i)
                                if sum(~inds_j)>1
                                    Decs(j,ii) = randsample(allowMode_i(~inds_j),1);
                                elseif sum(~inds_j)==1
                                    Decs(j,ii) = allowMode_i(~inds_j);
                                else
                                    Decs(j,ii) = randsample(allowMode_i,1);
                                end
                            end
%                             if i==N && length(unique(Decs(j,Decs(j,:)>0)))~=sum(Decs(j,:)>0)
%                                 error('Wrong Decisions');
%                             end
                        end
                    end
                    
                    Lower = repmat(obj.Global.upper_domain(1,:),size(Decs,1),1);
                    Upper = repmat(obj.Global.upper_domain(2,:),size(Decs,1),1);
                    Decs  = max(min(Decs,Upper),Lower);
                    
                case 'lower'
                    domain = obj.Global.lower_domain;
                    if ~iscell(domain)
                        Lower = repmat(domain(1,:),size(Decs,1),1);
                        Upper = repmat(domain(2,:),size(Decs,1),1);
                    else
                        Lower = domain{1};
                        Upper = domain{2};
                    end
                    Decs  = max(min(Decs,Upper),Lower);
            end
        end
        
        function [objs,adds] = CalObj(obj,Population,type)
            
            if isa(Population,'INDIVIDUAL')
                N = length(Population);
                s = transpose(Population.upper_decs);
                f = transpose(Population.lower_decs);
            elseif iscell(Population)
                s = Population{1};
                f = Population{2};
                N = size(s,2);
            end
            
            UE = obj.Parameter.UE;
            U = obj.Parameter.U;
            
            v_h = s(UE.N+1:2*UE.N,:);
            v_c = s(2*UE.N+1:end,:);
            s = s(1:UE.N,:);
            
            parameter = obj.Global.Parameter;
            
            R = zeros(UE.N,N);                  % uplink data rate
            R1 = R;
            
            for k=1:N
                %Resource allocation         
                Hsum = sum(U.H(s(:,k)==0,1));
                for i = 1:UE.N
                    if s(i,k)==0
                        R(i,k) = parameter.W*log2(1+UE.pt*U.H(i,1)/(parameter.wn+Hsum-U.H(i,1)));
                        R1(i,k) = inf;
                    elseif s(i,k)==i
                        R(i,k) = inf;
                        R1(i,k) = inf;
                    elseif s(i,k)>0
                        R(i,k) = parameter.W*log2(1+UE.pt*U.H(i,s(i,k)+1)/parameter.wn);
                        R1(i,k) = R(i,k);
                    elseif s(i,k)==-1
                        R(i,k) = inf;
                        R1(i,k) = inf;
                    end
                end
            end
            
%             T = repmat(U.F,1,N)./f(1:end-2,:) + (repmat(U.D,1,N)+repmat(U.B,1,N))./R;
            T = repmat(U.F,1,N)./f + (repmat(U.D,1,N)+repmat(U.B,1,N))./R;
            T(s==repmat(transpose(1:UE.N),1,N))=max(T(s==repmat(transpose(1:UE.N),1,N)),UE.Th);
            T(s~=repmat(transpose(1:UE.N),1,N))= T(s~=repmat(transpose(1:UE.N),1,N)) + UE.Th;
            obj.Parameter.T = T;
            
%             Ec = parameter.alpha.*f(1:end-2,:).^2.*repmat(U.F,1,N).*(s>0);
            Ec = parameter.alpha.*f.^2.*repmat(U.F,1,N).*(s>0);
            up = ((UE.pt.*repmat(U.D,1,N)+UE.pr.*repmat(U.B,1,N))./R).*(s~=repmat(transpose(1:UE.N),1,N));
            back = ((UE.pt.*repmat(U.B,1,N)+UE.pr.*repmat(U.D,1,N))./R1).*(s~=0 & s~=repmat(transpose(1:UE.N),1,N));
            
            pb = zeros(UE.N,N);
            EC = zeros(UE.N,N);
            ET_up = zeros(UE.N,N);
            ET_back = zeros(UE.N,N);
            for i = 1:N
                index = transpose(1:UE.N);
                ind = s(:,i)~= index;
                ET_up(index(ind),i) = up(ind,i);
                ind = s(:,i)== index;
                EC(s(ind,i),i) = Ec(ind,i);
                ind = (s(:,i)~= 0 & s(:,i)~= index);
                EC(s(ind,i),i) = Ec(ind,i);
                ET_back(s(ind,i),i) = back(ind,i);
                pb(:,i) = max([ET_up(:,i)+ET_back(:,i),EC(:,i)],[],2)./(UE.u*UE.Th*U.H(:,1));
            end
            
            adds = [T',EC',ET_up',ET_back',pb'];
            
            switch type
                case 'upper'
                    objs(:,1) = transpose(sum(T)/UE.N);
                    objs(:,2) = transpose(sum(f.*v_c*1e-9+pb.*v_h));
                case 'lower'
                    objs(:,1) = transpose(sum(Ec+up+back));
                    objs(:,2) = -transpose(sum(f.*v_c*1e-9+pb.*v_h-0.1*(pb*UE.Th+parameter.k*repmat(U.F,1,N).*(s==0))));
%                     objs(:,2) = -transpose(sum(f(1:end-2,:).*f(end-1,:)*1e-9+pb.*f(end,:)-0.1*(pb*UE.Th+parameter.k*repmat(U.D,1,N).*(s==0))));
           
            end
            
            
        end
        
        function cons = CalCon(obj,Population,type)
            
            if isa(Population,'INDIVIDUAL')
                s = transpose(Population.upper_decs);
                f = transpose(Population.lower_decs);
            elseif iscell(Population)
                s = Population{1};
                f = Population{2};
            end
            
            UE = obj.Parameter.UE;
            
            switch type
                case 'upper'
                    cons = transpose(obj.Parameter.T-UE.Tmax);
%                     cons(:,2) = transpose(sum(f.*(s(1:end-2,:)==0))-UE.f(1));
%                     cons(:,2) = transpose(max(f.*(s(1:end-2,:)>0)-UE.f(s(1:end-2,:)+1).*(s(1:end-2,:)>0)));
                case 'lower'
%                     cons(:,1) = transpose(max(obj.Parameter.T-UE.Tmax));
                    cons(:,1) = transpose(sum(f.*(s(1:UE.N,:)==0))-UE.f(1));
%                     cons(:,3) = transpose(max(f.*(s(1:end-2,:)>0)-UE.f(s(1:end-2,:)+1).*(s(1:end-2,:)>0)));
            end
        end
        
        function P = PF(obj)
            P=[];
% 			file = fullfile('Problems','DS','DS1.mat');
%             try
%                 load(file);
%             catch
%                 N = 1e3;
%                 r = obj.Parameter.('r');
%                 t=(2:1/(N-1):3)';
%                 P(:,1) = 1+r+(1+r)*cos(pi/2*t);
%                 P(:,2) = 1+r+(1+r)*sin(pi/2*t);
%                 save(file,'P');
%             end
        end
        
        function P = lower_PF(obj,y)
            P=[];
%             N = 1e6;
% %             r = obj.Parameter.('r');
%             t=(0:1/(N-1):y(:,1))';
%             P(:,1) = t.^2;
%             P(:,2) = (t-y(:,1)).^2;
        end
    end
    
end